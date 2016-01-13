#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Managing misfit windows
"""
import sys
import os.path
import re
#
import numpy as np
import scipy.signal as signal
#
import json
#
from obspy import UTCDateTime, read
from obspy.core.util.geodetics import gps2DistAzimuth, kilometer2degrees
from obspy.taup import TauPyModel
from obspy.imaging.beachball import Beach
#
import pyproj
#
from lanczos_interp1 import lanczos_interp1
#
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#
from taper import *


#====== utility functions
def is_equal(lst):
    return not lst or [lst[0]]*len(lst) == lst

#======
class Misfit(object):
    """Class managing all misfit windows

    data structure:

    Misfit: {
    |   'iteration':
    |   'events': {
    |   *   <event_id>: {
    |   |   |   'gcmt': {code, lat,lon,depth, ...}, 
    |   |   |   'relocate': {lat,lon,depth, ...}, 
    |   |   |   'stations': {
    |   |   |   *   <station_id(net.sta.loc)>: {
    |   |   |   |   |   'stat': {code, msg}, # code<0: problematic, code=0: OK
    |   |   |   |   |   'meta': {lat,lon,ele,...,channels:[{code,az,dip,...},...]},
    |   |   |   |   |   'filter': {type:, freqlim:},
    |   |   |   |   |   'taper': {type:, ratio:},
    |   |   |   |   |   'first_arrtime': , # used for cut noise window 
    |   |   |   |   |   'windows': {
    |   |   |   |   |   *   <window_id(cmp.pha)>: {
    |   |   |   |   |   |   |   component, azimuth, dip, starttime, endtime, weight,
    |   |   |   |   |   |   |   'stat': {code, msg}, #code<0: problematic, code=0: not measured, code=1: measured
    |   |   |   |   |   |   |   'phase': {name, ttime, takeoff_angle, ray_param},
    |   |   |   |   |   |   |   'quality': {A_obs, A_noise, SNR, A_syn},
    |   |   |   |   |   |   |   'misfit': {CC0, CCmax, CC_time_shift,
    |   |   |   |   |   |   |              AR0, ARmax } },
    |   |   |   |   |   *   <window_id>: { },
    |   |   |   |   |   |   ...
    |   |   |   *   <station_id>: { ... },
    |   |   |   |   ...
    |   *   <event_id>: { ... },
    |   |   ...

    NOTE:
        1. 1D Earth model: ak135
    """

    def __init__(self):
        """Misfit dict
        """
        self.data = {'events':{}}


    def save(self, filename='misfit.json'):
        """Save data
        """
        with open(filename, 'w') as fp:
            json.dump(self.data, fp, indent=2)
            

    def load(self, filename='misfit.json'):
        """Load data
        """
        with open(filename, 'r') as fp:
            self.data = json.load(fp)


    def setup_event_from_CMTSOLUTION(self, cmt_file, update=False):
        """cmt_file (str): CMTSOLUTION format file
        """
        with open(cmt_file, 'r') as f:
            cmt = [ x.split() for x in f.readlines() 
                    if not(x.startswith('#')) ]
        year  = cmt[0][1]
        month  = cmt[0][2]
        day = cmt[0][3]
        hour  = cmt[0][4]
        minute = cmt[0][5]
        second = cmt[0][6]

        event_id = cmt[1][2]
        time_shift = float(cmt[2][2])
        half_duration = float(cmt[3][2])
        lat = float(cmt[4][1])
        lon = float(cmt[5][1])
        depth = float(cmt[6][1])

        # centroid time 
        isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(
                year, month, day, hour, minute, second)
        centroid_time = UTCDateTime(isotime) + time_shift

        # moment tensor
        # basis: (r,theta,phi) corresponds to (up,south,east)
        Mrr = float(cmt[7][1])
        Mtt = float(cmt[8][1])
        Mpp = float(cmt[9][1])
        Mrt = float(cmt[10][1])
        Mrp = float(cmt[11][1])
        Mtp = float(cmt[12][1])
        M = [[Mrr, Mrt, Mrp], [Mrt, Mtt, Mtp], [Mrp, Mtp, Mpp]]

        # add event
        events = self.data['events']
        gcmt = {'code': event_id,
                'centroid_time':str(centroid_time),
                'half_duration':half_duration,
                'latitude':lat, 'longitude':lon, 'depth':depth, 
                'moment_tensor':M }
        if event_id not in events:
            events[event_id] = {
                    'gcmt': gcmt,
                    'stations': {},
                    'stat': {'code': 0, 'msg': ""}}
        elif update:
            events[event_id]['gcmt'].update(gcmt)
            events[event_id]['stat']['msg'] = "updated"
            print "[WARNING] %s: update event info" % (event_id)
        else:
            print "[WARNING] %s: already exists, skip." % (event_id)


    def setup_stations_from_metafile(self, station_file,
            event_id_list=None, band_code=None, three_channels=True, 
            update=False):
        """ Setup station metadata. 

            station_file (str): 
                FDSN-station text format file at channel level
            event_id_list (list): 
                list of event ID's to which stations are added [default: None]
                default to all events.
            band_code (str):
                instrument/band code [default: None]
            three_channels (bool):
                check for completeness of 3 channels [default: False]

            Note: only use stations which have the same lat/lon/ele/depth 
                in all the available channels.
        """
        # read station file
        with open(station_file, 'r') as f:
            lines = [x.replace('\n','').split('|')  \
                    for x in f.readlines() if not(x.startswith('#'))]
        
        # get all station metadata
        stations_all = {}
        for x in lines:
            net_sta_loc = (x[0], x[1], x[2])
            date1 = [ int(a) for a in re.sub("\D", " ", x[15]).split() ]
            date2 = [ int(a) for a in re.sub("\D", " ", x[16]).split() ]
            t1 = UTCDateTime(date1[0], date1[1], date1[2]) \
                    + 60.0*(60.0*date1[3] + date1[4]) + date1[5]
            t2 = UTCDateTime(date2[0], date2[1], date2[2]) \
                    + 60.0*(60.0*date2[3] + date2[4]) + date2[5]
            channel = {'code':        x[3],
                       'latitude':    float(x[4]),
                       'longitude':   float(x[5]),
                       'elevation':   float(x[6]),
                       'depth':       float(x[7]),
                       'azimuth':     float(x[8]),
                       'dip':         float(x[9]),
                       'starttime':   t1.isoformat(),
                       'endtime':     t2.isoformat() }
            if net_sta_loc not in stations_all:
                stations_all[net_sta_loc] = []
            stations_all[net_sta_loc].append(channel)

        # empty/Null event_id_list defaults to all events
        events = self.data['events']
        if not event_id_list:
            event_id_list = [ x for x in events ]

        # add stations for each event
        for event_id in event_id_list:
            # check if event info is set 
            if event_id not in events:
                print "[WARNING] %s does not exist, SKIP" % (event_id)
                continue

            event = events[event_id]
            stations = event['stations']
            gcmt = event['gcmt']

            # station active time is set to centroid time
            active_time = UTCDateTime(gcmt['centroid_time'])

            for net_sta_loc in stations_all:
                # station_id: net.sta.loc
                station_id = '.'.join(net_sta_loc)

                # skip existing stations if not update
                if (station_id in stations) and (not update):
                    print "[WARNING] %s already exists, SKIP." \
                            % (station_id)
                    continue

                # select channels which are active at the specified time 
                channels = [ x for x in stations_all[net_sta_loc] 
                        if UTCDateTime(x['starttime']) < active_time 
                        and UTCDateTime(x['endtime']) > active_time ]
                # select band code
                if band_code:
                    n = len(band_code)
                    channels = [ x for x in channels 
                            if x['code'][0:n]==band_code ]
                # check same locations for all selected channels 
                lats = [ x['latitude'] for x in channels ] 
                lons = [ x['longitude'] for x in channels ] 
                eles = [ x['elevation'] for x in channels ] 
                deps = [ x['depth'] for x in channels ] 
                if lats.count(lats[0])!=len(lats) or \
                        lons.count(lons[0])!=len(lons) or \
                        eles.count(eles[0])!=len(eles) or \
                        deps.count(deps[0])!=len(deps):
                    print "[WARNING] %s: " \
                          "channels do NOT have the same coordinates, SKIP." \
                          % (station_id)
                    continue
                # check completeness for 3-components
                if three_channels:
                    if len(channels) != 3:
                        print '[WARNING] %s: not exactly 3 components found, '\
                              'SKIP.' % (station_id)
                        continue
                    # check channel orientations
                    Z_comp = [ (x['code'], x['azimuth'], x['dip'])
                            for x in channels if x['code'][2] == 'Z']
                    if len(Z_comp) != 1 or abs(Z_comp[0][2]) != 90.0: 
                        print '[WARNING] %s: problematic Z channel, SKIP' \
                                % (station_id)
                        print '          channel: ', H_comp
                        continue
                    H_comp = [ (x['code'], x['azimuth'], x['dip']) \
                            for x in channels if x['code'][2] != 'Z']
                    if len(H_comp) != 2 or \
                            abs(H_comp[0][2]) != 0.0 or \
                            abs(H_comp[1][2]) != 0.0 or \
                            abs(np.cos(np.deg2rad(
                                H_comp[0][1] - H_comp[1][1]))) > 0.1: 
                        print '[WARNING] %s: problematic horizontal channels, SKIP'\
                                % (station_id)
                        print '          channel: ', H_comp
                        continue

                # geodetic
                dist, az, baz = gps2DistAzimuth(gcmt['latitude'], gcmt['longitude'],
                        channels[0]['latitude'], channels[0]['longitude'])
                dist_degree = kilometer2degrees(dist/1000.0)

                # form station metadata 
                meta = {'latitude': channels[0]['latitude'], 
                        'longitude': channels[0]['longitude'],
                        'elevation': channels[0]['elevation'],
                        'depth': channels[0]['depth'],
                        'azimuth': az,
                        'back_azimuth': baz,
                        'dist_degree': dist_degree,
                        'channels': channels }

                # add station info
                if station_id not in stations:
                    stations[station_id] = {
                            'meta': meta,
                            'windows': {},
                            'stat': {'code': 0, 'msg': ""} }
                elif update:
                    stations[station_id]['meta'].update(meta)
                    stations[station_id]['stat']['msg'] = "updated"

            #end for net_sta_loc in stations_all:
        #end for event_id in event_id_list:


    def setup_windows(self,
            window_list=[('F','p,P',[-30,50]), ('F','s,S',[-40,70])],
            filter_param=('butter', 3, [0.012, 0.08]),
            taper_param=('cosine', 0.1),
            event_id_list=None, station_id_list=None):
        """window_list: [ (component, phases, [begin, end]), ...]
           filter_param: (type, freqlims)
        """
        # check event/station_id_list
        events = self.data['events']
        if not event_id_list:
            event_id_list = [ x for x in events ]
        if not station_id_list:
            loop_all_stations = True
        else:
            loop_all_stations = False

        # initiate taup
        taup_model = TauPyModel(model="ak135")

        # filter/taper parameters
        filter_dict = {'type': filter_param[0], 
                'order': filter_param[1], 'freqlim': filter_param[2]}
        taper_dict = {'type': taper_param[0], 'ratio': taper_param[1]}
        if not 0.0 < taper_param[1] < 0.5:
            raise ValueError("taper ratio must lie between 0 and 0.5.")

        # loop each event
        for event_id in event_id_list:
            if event_id not in events:
                print "[WARNING] %s does NOT exist. SKIP" \
                        % (event_id)
                continue

            event = events[event_id]
            gcmt = event['gcmt']
            centroid_time = UTCDateTime(gcmt['centroid_time'])
            stations = event['stations']
            if loop_all_stations:
                station_id_list = [ x for x in stations ]

            # loop each station
            for station_id in station_id_list:
                if station_id not in stations:
                    print "[WARNING] %s:%s does NOT exist. SKIP" \
                            % (event_id, station_id)
                    continue

                station = stations[station_id]
                meta = station['meta']
                baz = meta['back_azimuth']
                dist_degree = meta['dist_degree']

                # setup filter/taper parameters
                station['filter'] = filter_dict
                station['taper'] = taper_dict

                # setup first arrival time 
                arrivals = taup_model.get_travel_times(
                        source_depth_in_km=gcmt['depth'],
                        distance_in_degree=dist_degree,
                        phase_list=['ttp'])
                if arrivals:
                    arr = arrivals[0]
                    first_arrtime = str(centroid_time + arr.time)
                else:
                    print "[WARNING] %s:%s ttp not found, skip." \
                            % (event_id, station_id)
                    continue
                station['first_arrtime'] = first_arrtime 

                # setup singal windows
                windows = station['windows']
                for win in window_list:
                    comp = win[0]
                    phase = win[1]
                    signal_begin = float(win[2][0])
                    signal_end = float(win[2][1])
                    window_id = "%s.%s" % (comp, phase)

                    if comp == 'Z': # vertcal component
                        cmpaz = 0.0 
                        cmpdip = -90.0
                    elif comp == 'R': # radial component
                        cmpaz = (baz + 180.0)%360.0
                        cmpdip = 0.0
                    elif comp == 'T': # tangential component (TRZ: right-hand convention)
                        cmpaz = (baz - 90.0)%360.0
                        cmpdip = 0.0
                    elif comp == 'H': # horizontal vector
                        cmpaz = float('nan')
                        cmpdip = 0.0
                    elif comp == 'F': # full vector
                        cmpaz = float('nan')
                        cmpdip = float('nan')
                    else:
                        print "[WARNING] %s: unrecognized component, SKIP." \
                                % (comp)
                        continue
                    # phase arrivals predicted from reference Earth model
                    arrivals = taup_model.get_travel_times(
                            source_depth_in_km=gcmt['depth'],
                            distance_in_degree=dist_degree,
                            phase_list=phase.split(','))
                    if not arrivals:
                        print "[WARNING] %s:%s:%s phase not found, skip " \
                                % (event_id, station_id, phase)
                        continue
                    # use first arrival as the window reference time
                    arr = arrivals[0]
                    signal_starttime = str(centroid_time + arr.time + signal_begin)
                    signal_endtime = str(centroid_time + arr.time + signal_end)
                    # add window
                    window = {
                        'component': comp,
                        'azimuth': cmpaz,
                        'dip': cmpdip,
                        'starttime': signal_starttime,
                        'endtime': signal_endtime,
                        'phase': {'name': arr.name, 'ttime': arr.time,
                                  'takeoff_angle': arr.takeoff_angle,
                                  'ray_param': arr.ray_param},
                        'quality': {},
                        'misfit': {},
                        'stat': {'code': 0, 'msg': ""} }
                    windows[window_id] = window

                #end for winpar in winpar_list:
            #end for station_id in station_id_list:
        #end for event_id in event_id_list:


    def read_seismograms(self, event_id, station_id, 
            obs_dir='obs', syn_dir='syn', 
            syn_band_code='MX', syn_suffix='.sem.sac'):
        """ Get 3-component seismograms recorded at one station for one event.
            Return:
                syn_st, obs_st
        """
        syn_orientation_codes = ['E', 'N', 'Z']
        events = self.data['events']
        event = events[event_id]
        stations = event['stations']
        station = stations[station_id]
        meta = station['meta']
        channels = meta['channels']
        first_arrtime = UTCDateTime(station['first_arrtime'])
        windows = station['windows']
        window_endtime = [UTCDateTime(windows[x]['endtime']) for x in windows]
        max_window_endtime = max(window_endtime)
 
        #------ get file paths of obs, syn seismograms
        obs_files = [ '{:s}/{:s}.{:s}'.format(
            obs_dir, station_id, x['code']) for x in channels ]
        syn_files = [ '{:s}/{:s}.{:2s}{:1s}{:s}'.format(
            syn_dir, station_id, syn_band_code, x, syn_suffix)
            for x in syn_orientation_codes ]

        #------ read in obs, syn seismograms
        obs_st  = read(obs_files[0])
        obs_st += read(obs_files[1])
        obs_st += read(obs_files[2])
        syn_st  = read(syn_files[0])
        syn_st += read(syn_files[1])
        syn_st += read(syn_files[2])

        #------ get time samples of syn seismograms
        if not is_equal( [ (tr.stats.starttime, tr.stats.delta, tr.stats.npts) \
                for tr in syn_st ] ):
            raise Exception('%s:%s: not equal time samples in'\
                    ' synthetic seismograms.' % (event_id, station_id))
        tr = syn_st[0]
        syn_starttime = tr.stats.starttime
        syn_endtime = tr.stats.endtime
        syn_delta = tr.stats.delta
        syn_npts = tr.stats.npts
        syn_times = syn_delta*np.arange(syn_npts)

        #------ interpolate obs into the same time samples of syn
        obs_ENZ = np.zeros((3, syn_npts))
        #taper_len = 50.0 # FIXME this value is hard wired now
        syn_nyq = 0.5/syn_delta
        for i in range(3):
            tr = obs_st[i]
            obs_npts = tr.stats.npts
            obs_delta = tr.stats.delta
            obs_starttime = tr.stats.starttime
            obs_endtime = tr.stats.endtime
            # check if obs record is long enough (60.0 sec more)
            if (first_arrtime - obs_starttime) < 60.0 or \
               (obs_endtime - max_window_endtime) < 60.0:
                raise Exception('%s:%s: obs does not start 60s before first arrtime ' \
                        'or end 60s after the last window endtime' \
                        % (event_id, station_id))
            # lowpass to frequency band in sampling rate of syn
            tr.detrend(type='linear')
            #tr.taper(max_percentage=0.1, max_length=taper_len)
            tr.filter('lowpassCheby2', freq=0.9*syn_nyq, maxorder=12)
            # interpolation
            obs_ENZ[i,:] = lanczos_interp1(tr.data, obs_delta,
                    syn_times+(syn_starttime-obs_starttime), na=20)

        #------ rotate obs to ENZ
        # projection matrix: obs = proj * ENZ => ZNE = inv(proj) * obs
        proj_matrix = np.zeros((3, 3))
        for i in range(3):
            channel = channels[i]
            sin_az = np.sin(np.deg2rad(channel['azimuth']))
            cos_az = np.cos(np.deg2rad(channel['azimuth']))
            sin_dip = np.sin(np.deg2rad(channel['dip']))
            cos_dip = np.cos(np.deg2rad(channel['dip']))
            # column vector = obs channel polarization 
            proj_matrix[i,0] = cos_dip*sin_az # proj to E
            proj_matrix[i,1] = cos_dip*cos_az # proj to N
            proj_matrix[i,2] = -sin_dip       # proj to Z
        # inverse projection matrix: ENZ = inv(proj) * obs
        inv_proj = np.linalg.inv(proj_matrix)
        obs_ENZ = np.dot(inv_proj, obs_ENZ)

        #------ form obs traces
        obs_st = syn_st.copy()
        for i in range(3):
            obs_st[i].data = obs_ENZ[i,:]

        return syn_st, obs_st


    def measure_windows_for_one_station(self, event_id, station_id, 
            obs_dir='obs', syn_dir='syn',
            syn_band_code='MX', syn_suffix='.sem.sac',
            plot=True, use_STF=False,
            cc_delta=0.01, output_adj=False, adj_dir='adj',
            output_hess=False, hess_dir='hess',
            adj_window_id_list=['F.p,P', 'F.s,S'],
            weight_param={'SNR':[10,15], 'CCmax':[0.6,0.8], 'CC0':[0.5,0.7]}):
        """
            output_hess, hess_dir: 
               compute adjoint sources for approximating diagonal hess matrix.
        """
        #
        events = self.data['events']
        event = events[event_id]
        cmt = event['gcmt']
        half_duration = cmt['half_duration']
        centroid_time = UTCDateTime(cmt['centroid_time'])
        stations = event['stations']
        station = stations[station_id]
        windows = station['windows']
        syn_orientation_codes = ['E', 'N', 'Z']
        first_arrtime = UTCDateTime(station['first_arrtime'])

        #------ read obs/syn seismograms
        # syn: u, obs: d
        syn_st, obs_st = self.read_seismograms(
                event_id, station_id,
                obs_dir=obs_dir, syn_dir=syn_dir,
                syn_band_code=syn_band_code, syn_suffix=syn_suffix)

        # time samples of syn 
        tr = syn_st[0]
        syn_starttime = tr.stats.starttime
        syn_delta = tr.stats.delta
        syn_sampling_rate = tr.stats.sampling_rate
        syn_npts = tr.stats.npts
        syn_times = syn_delta*np.arange(syn_npts)
        syn_nyq = 0.5/syn_delta

        # syn/obs data arrays
        syn_ENZ = np.zeros((3, syn_npts))
        obs_ENZ = np.zeros((3, syn_npts))
        for i in range(3):
            syn_ENZ[i,:] = syn_st[i].data
            obs_ENZ[i,:] = obs_st[i].data

        # filter parameters
        filter_param = station['filter']
        filter_type = filter_param['type']
        filter_order = filter_param['order']
        filter_freqlim = filter_param['freqlim']

        # window taper parameters
        taper_param = station['taper']
        taper_type = taper_param['type']
        taper_ratio = taper_param['ratio']

        #------ filter spectrum
        # F = filter of data measurment
        filter_b, filter_a = signal.butter(filter_order,
                np.array(filter_freqlim)/syn_nyq, btype='band')
        if use_STF:
            # S: source spectrum (isosceles triangle)
            f = np.fft.rfftfreq(syn_npts, d=syn_delta)
            F_src = np.sinc(f * half_duration)**2

        # filter syn/obs
        # obs = F * d
        obs_ENZ[:,:] = signal.filtfilt(filter_b, filter_a, obs_ENZ)
        # syn = F * S * u
        syn_ENZ[:,:] = signal.filtfilt(filter_b, filter_a, syn_ENZ)
        if use_STF:
            syn_ENZ[:,:] = np.fft.irfft(F_src*np.fft.rfft(syn_ENZ), syn_npts)

        # noise: use signals 40s before first arrival time on obs
        #FIXME: better choice of the time length before first arrival? 
        t0 = (first_arrtime - syn_starttime) - 40.0
        noise_idx = syn_times < t0
        #t = syn_times[noise_idx]
        #b = t[0]
        #e = t[-1]
        #taper_width = (e-b) * 0.1
        #win_c = [b, b+taper_width, e-taper_width, e]
        #taper = cosine_taper(t, win_c)
        noise_ENZ = obs_ENZ[:,noise_idx]

        #------ process windows
        if output_adj:
            adj_ENZ = np.zeros((3, syn_npts))
        if output_hess:
            hess_ENZ = np.zeros((3, syn_npts))
        taper = np.zeros(syn_npts)
        taper = np.zeros(syn_npts)
        proj_matrix = np.zeros((3,3))
        cc = np.zeros(2*syn_npts-1)
        for window_id in windows:
            # window parameters
            window = windows[window_id]
            comp = window['component']
            cmpaz = window['azimuth']
            cmpdip = window['dip']
            window_starttime = UTCDateTime(window['starttime'])
            window_endtime = UTCDateTime(window['endtime'])
            window_len = window_endtime - window_starttime
            # window taper
            win_b = window_starttime - syn_starttime
            win_e = window_endtime - syn_starttime
            taper_width = window_len * min(taper_ratio, 0.5)
            win_c = [win_b, win_b+taper_width, win_e-taper_width, win_e]
            taper = cosine_taper(syn_times, win_c)
            # projection matrix
            proj_matrix[:,:] = 0.0 #reset to zero
            if comp in ['Z', 'R', 'T']:
                sin_az = np.sin(np.deg2rad(cmpaz))
                cos_az = np.cos(np.deg2rad(cmpaz))
                sin_dip = np.sin(np.deg2rad(cmpdip))
                cos_dip = np.cos(np.deg2rad(cmpdip))
                n = np.array([ [cos_dip * sin_az], # cos(E, comp)
                               [cos_dip * cos_az], # N, comp
                               [-sin_dip] ])       # Z, comp
                proj_matrix = np.dot(n, n.transpose())
            elif comp == 'H': # horizontal vector 2d
                proj_matrix[0,0] = 1.0 # E
                proj_matrix[1,1] = 1.0 # N
                proj_matrix[2,2] = 0.0 # Z
            elif comp == 'F': # full 3d vector
                proj_matrix[0,0] = 1.0
                proj_matrix[1,1] = 1.0
                proj_matrix[2,2] = 1.0
            else:
                print '[WARNING] %s:%s:%s unrecognized component code, SKIP' \
                        % (event_id, station_id, window_id)
                continue

            # apply window taper and projection
            # obs = w * F * d
            obs_ENZ_win = np.dot(proj_matrix, obs_ENZ) * taper
            # obs = w * F * S * u
            syn_ENZ_win = np.dot(proj_matrix, syn_ENZ) * taper
            # noise
            noise_ENZ_win = np.dot(proj_matrix, noise_ENZ)

            # measure SNR
            A_syn = np.sqrt(np.max(np.sum(syn_ENZ_win**2, axis=0)))
            A_obs = np.sqrt(np.max(np.sum(obs_ENZ_win**2, axis=0)))
            A_noise =  np.sqrt(np.max(np.sum(noise_ENZ_win**2, axis=0)))
            if A_obs == 0: # data file may not containe the time window
                print '[WARNING] %s:%s:%s empty obs trace, SKIP.' \
                        % (event_id, station_id, window_id)
                continue
            if A_noise == 0: # could occure when the data begin time is too close to the first arrival
                print '[WARNING] %s:%s:%s empty noise trace, SKIP.' \
                        % (event_id, station_id, window_id)
                continue
            snr = 20.0 * np.log10(A_obs/A_noise)
 
            # normalized cross-correlation
            obs_norm2 = np.sum(obs_ENZ_win**2)
            obs_norm = np.sqrt(obs_norm2)
            syn_norm = np.sqrt(np.sum(syn_ENZ_win**2))
            cc[:] = 0.0
            for i in range(3):
                # NOTE the order (obs,syn) is important. The positive time on 
                #   CC means shifting syn in the positive time direction
                cc += signal.fftconvolve(
                        obs_ENZ_win[i,:], syn_ENZ_win[i,::-1], 'full')
            cc /= obs_norm * syn_norm
            # zero-lag cc coeff.
            CC0 = cc[syn_npts]
            AR0 = CC0 * syn_norm / obs_norm # amplitude ratio syn/obs 

            # cross-correlation time shift
            # tshift>0: synthetic is shifted along the positive time direction
            CC_shift_range = window_len/2.0 #TODO: more reasonable choice?
            ncc = int(CC_shift_range / cc_delta)
            cc_times = np.arange(-ncc,ncc+1) * cc_delta
            # interpolate cc to finer time samples
            if syn_delta < cc_delta:
                print '[WARNING] syn_delta(%f) < cc_time_step(%f)' \
                        % (syn_delta, cc_delta)
            ti = cc_times + (syn_npts-1)*syn_delta  # -(npts-1)*dt: begin time in cc
            cci = lanczos_interp1(cc, syn_delta, ti, na=20)
            # time shift at the maximum correlation
            imax = np.argmax(cci)
            CC_time_shift = cc_times[imax]
            CCmax = cci[imax]
            ARmax = CCmax * syn_norm / obs_norm # amplitude ratio: syn/obs

            # window weighting based on SNR and misfit
            weight = cosine_taper(snr, weight_param['SNR']) \
                    * cosine_taper(CCmax, weight_param['CCmax']) \
                    * cosine_taper(CC0, weight_param['CC0'])

            # adjoint source (objective functional: zero-lag cc coef.)
            # adj = conj(F * S) * w * [ w * F * d - A * w * F * S * u] / N, 
            #   where A = CC0(un-normalized) / norm(syn)**2, N = norm(obs)*norm(syn)
            if output_adj:
                A0 = CC0 * obs_norm / syn_norm # amplitude raito: obs/syn
                adj =  taper * (obs_ENZ_win - A0*syn_ENZ_win) / obs_norm / syn_norm
                adj_ENZ_win = signal.filtfilt(filter_b, filter_a, adj)
                A_adj = np.sqrt(np.max(np.sum(adj_ENZ_win**2, axis=0)))
                if use_STF:
                    adj_ENZ_win = np.fft.irfft(np.conjugate(F_src) * 
                            np.fft.rfft(adj_ENZ_win), syn_npts)
                if window_id in adj_window_id_list: 
                    adj_ENZ += weight * adj_ENZ_win

            # adjoint source to approximate hessian diagonal
            # hess = conj(F * S) * w * [w * F * S * u] / N, 
            #   where N = norm(syn)**2
            if output_hess:
                hess =  taper * syn_ENZ_win / (syn_norm**2)
                hess_ENZ_win = signal.filtfilt(filter_b, filter_a, hess)
                if use_STF:
                    hess_ENZ_win = np.fft.irfft(np.conjugate(F_src) * 
                            np.fft.rfft(hess_ENZ_win), syn_npts)
                if window_id in adj_window_id_list: 
                    hess_ENZ += weight * hess_ENZ_win

            # record results
            quality_dict = {
                    'A_obs': A_obs, 'A_syn': A_syn, 'A_noise': A_noise,
                    'SNR': snr}
            misfit_dict = {
                    'CC0': CC0, 'CCmax': CCmax,
                    'AR0': AR0, 'ARmax': ARmax,
                    'CC_time_shift': CC_time_shift }
            window['quality'] = quality_dict
            window['misfit'] = misfit_dict
            window['weight'] = weight
            window['stat'] = {'code': 1, 'msg': "measured"}

            # plot measure window and results 
            if plot:
                t = syn_times
                for i in range(3):
                    plt.subplot(411+i)
                    if i == 0:
                        plt.title('%s dt %.2f CCmax %.3f ARmax %.3f CC0 %.3f '
                                'AR0 %.3f \nAobs %g Anoise %g SNR %.1f weight %.3f'
                                % (station_id, CC_time_shift, CCmax, ARmax, 
                                    CC0, AR0, A_obs, A_noise, snr, weight) )
                    plt.plot(t, obs_ENZ[i,:]/A_obs, 'k', linewidth=0.2)
                    plt.plot(t, syn_ENZ[i,:]/A_syn, 'r', linewidth=0.2)
                    plt.plot(t[noise_idx], noise_ENZ_win[i,:]/A_obs, 'b', linewidth=1.0)
                    idx = (win_b <= syn_times) & (syn_times <= win_e)
                    plt.plot(t[idx], obs_ENZ_win[i,idx]/A_obs, 'k', linewidth=1.0)
                    plt.plot(t[idx], syn_ENZ_win[i,idx]/A_obs * A0, 'r', linewidth=1.0)
                    if output_adj:
                        plt.plot(t, adj_ENZ_win[i,:]/A_adj, 'c', linewidth=1.0)
                    plt.ylim((-1.5, 1.5))
                    plt.xlim((min(t), max(t)))
                    plt.ylabel(syn_orientation_codes[i])
                plt.subplot(414)
                plt.plot(cc_times, cci, 'k-')
                plt.xlim((min(cc_times), max(cc_times)))
                plt.ylabel(window_id)
                plt.show()
        #end for window_id in windows:

        # output adjoint source for the current station
        if output_adj:
            for i in range(3):
                tr = syn_st[i]
                tr.data = adj_ENZ[i,:]
                out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
                        adj_dir, station_id, syn_band_code, 
                        syn_orientation_codes[i])
                # sac format 
                tr.write(out_file + '.adj.sac', 'sac')
                # ascii format 
                # time is relative to event origin time
                shift = tr.stats.sac.b - tr.stats.sac.o
                with open(out_file+'.adj','w') as fp:
                    for j in range(syn_npts):
                        fp.write("{:16.9e}  {:16.9e}\n".format(
                            syn_times[j]+shift, adj_ENZ[i,j]))

        # output hess adjoint source for the current station
        if output_hess:
            for i in range(3):
                tr = syn_st[i]
                tr.data = hess_ENZ[i,:]
                out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
                        hess_dir, station_id, syn_band_code, 
                        syn_orientation_codes[i])
                # sac format 
                tr.write(out_file + '.adj.sac', 'sac')
                # ascii format 
                # time is relative to event origin time
                shift = tr.stats.sac.b - tr.stats.sac.o
                with open(out_file+'.adj','w') as fp:
                    for j in range(syn_npts):
                        fp.write("{:16.9e}  {:16.9e}\n".format(
                            syn_times[j]+shift, hess_ENZ[i,j]))


    def measure_windows_for_one_event(self, event_id,
            obs_dir='obs', syn_dir='syn',
            syn_band_code='MX', syn_suffix='.sem.sac',
            plot=False, use_STF=False,
            cc_delta=0.01, output_adj=False, adj_dir='adj', 
            output_hess=False, hess_dir='adj',
            adj_window_id_list=['F.p,P', 'F.s,S'],
            weight_param={'SNR':[10,15], 'CCmax':[0.6,0.8], 'CC0':[0.5,0.7]}):
        """measure misfit on time windoes for one event
            cc_delta: sampling interval for cross-correlation between obs and syn.
            weight_param: one-sided cosine taper [stop, pass]
            use_STF:
                flag to convolve source time function to synthetics. STF is 
                specifed as isosceles triangle with half_duration.
        """
        # check inputs
        events = self.data['events']
        if event_id not in events:
            print "[ERROR] %s does not exist."  % (event_id)
            sys.exit()
        event = events[event_id]

        # process each station
        stations = event['stations']
        for station_id in stations:
            station = stations[station_id]
            try:
                self.measure_windows_for_one_station(event_id, station_id, 
                    obs_dir=obs_dir, syn_dir=syn_dir, 
                    syn_band_code=syn_band_code, syn_suffix=syn_suffix, 
                    plot=plot, use_STF=use_STF,
                    cc_delta=cc_delta, output_adj=output_adj, adj_dir=adj_dir,
                    output_hess=output_hess, hess_dir=hess_dir,
                    adj_window_id_list=adj_window_id_list, 
                    weight_param=weight_param)
            except Exception as e:
                print "[WARNING] %s" % (str(e))
                station['stat']['code'] = -1
                station['stat']['msg'] = str(e)
                continue
            # log status
            station['stat']['code'] = 1
        #end for station_id in station_id_list:


    def relocate_1d(self, event_id,
            window_id_list=['F.p,P', 'F.s,S'], 
            fix_depth=False,
            out_cmt_file=None):
        """relocate event using ray path in reference earth model
        """
        # check inputs
        events = self.data['events']
        if event_id not in events:
            print "[ERROR] %s does NOT exist. Exit" % (event_id)
            sys.exit()

        # select windows
        sta_win_id_list = []
        event = events[event_id]
        stations = event['stations']
        for station_id in stations:

            station = stations[station_id]
            if station['stat']['code'] < 0:
                continue

            windows = station['windows']
            for window_id in window_id_list:
                if window_id not in windows:
                    continue

                window = windows[window_id]
                if window['stat']['code'] != 1:
                    continue 

                misfit = window['misfit']
                #if window['quality']['SNR'] < min_SNR or \
                #        misfit['CC0'] < min_CC0 or \
                #        misfit['CCmax'] < min_CCmax or\
                #        abs(misfit['CC_time_shift']) > max_CC_time_shift:
                #    continue

                sta_win_id = (station_id, window_id)
                sta_win_id_list.append(sta_win_id)

        # create sensitivity matrix G in local NED coordinate
        # G * dm  = dt_cc
        # G: [[-px_1, -py_1, -pz_1, 1.0], # ray1
        #     [-px_2, -py_2, -pz_2, 1.0], # ray2
        #     ...]
        # dm: [dNorth(km), dEast, dDepth, dT(sec)]
        # dt_cc: [dt1, dt2, ...]
        n = len(sta_win_id_list)
        G = np.zeros((n, 4))
        dt_cc = np.zeros(n)
        R_Earth_km = 6371.0
        gcmt = event['gcmt']
        evdp = gcmt['depth']
        for i in range(n):
            sta_win_id = sta_win_id_list[i]
            station_id = sta_win_id[0]
            window_id = sta_win_id[1]
    
            station = stations[station_id]
            meta = station['meta']
            window = station['windows'][window_id]
            phase = window['phase']
            misfit = window['misfit']
            weight = window['weight']

            azimuth = np.deg2rad(meta['azimuth'])
            takeoff_angle = phase['takeoff_angle']
            takeoff = np.deg2rad(takeoff_angle + 180.0*(takeoff_angle<0))
            ray_param = phase['ray_param']
            slowness = ray_param / (R_Earth_km - evdp) #unit: s/km
            # local coordinate: NED
            pd = np.cos(takeoff) * slowness
            pn = np.cos(azimuth) * np.sin(takeoff) * slowness
            pe = np.sin(azimuth) * np.sin(takeoff) * slowness
            # create sensitivity matrix 
            G[i,:] = weight * np.array([-pn, -pe, -pd, 1.0]) # -p: from receiver to source
            dt_cc[i] = weight * misfit['CC_time_shift']

        #linearized inversion (can be extended to second order using dynamic ray-tracing)
        if fix_depth: 
            G[:, 2] = 0.0
        dm, residual, rank, sigval = np.linalg.lstsq(G, dt_cc)

        # convert dm from NED to ECEF coordinate
        evla = gcmt['latitude']
        evlo = gcmt['longitude']

        slat = np.sin(np.deg2rad(evla))
        clat = np.cos(np.deg2rad(evla))
        slon = np.sin(np.deg2rad(evlo))
        clon = np.cos(np.deg2rad(evlo))

        N = [-slat*clon, -slat*slon, clat]
        E = [-slon, clon, 0.0]
        D = [-clat*clon, -clat*slon, -slat]

        ev_dx = N[0]*dm[0] + E[0]*dm[1] + D[0]*dm[2]
        ev_dy = N[1]*dm[0] + E[1]*dm[1] + D[1]*dm[2]
        ev_dz = N[2]*dm[0] + E[2]*dm[1] + D[2]*dm[2]
        ev_dt = dm[3]

        # initialize pyproj objects
        geod = pyproj.Geod(ellps='WGS84')
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

        # old location in ECEF (meters)
        evx, evy, evz = pyproj.transform(lla, ecef, evlo, evla, -1000.0*evdp)

        # new location in ECEF (meters)
        evx1 = evx + ev_dx*1000.0
        evy1 = evy + ev_dy*1000.0
        evz1 = evz + ev_dz*1000.0
        # in LLA
        evlo1, evla1, evalt1 = pyproj.transform(ecef, lla, evx1, evy1, evz1)
        evdp1 = -evalt1/1000.0

        # residuals 
        # linearized modelling
        dt_syn = G.dot(dm)
        dt_res = dt_cc - dt_syn

        # make results
        new_centroid_time = UTCDateTime(gcmt['centroid_time']) + ev_dt
        reloc_dict = {
                'window_id_list': window_id_list,
                'singular_value': sigval.tolist(),
                'dm': {'dNorth':dm[0], 'dEast':dm[1], 'dDepth':dm[2], 
                    'dT':dm[3]},
                'latitude':evla1, 
                'longitude':evlo1, 
                'depth':evdp1,
                'centroid_time': str(new_centroid_time),
                'data': {'num':n, 'mean':np.mean(dt_cc), 'std':np.std(dt_cc)},
                'residual': {'mean':np.mean(dt_res), 'std':np.std(dt_res)} }

        event['relocate'] = reloc_dict

        # make new CMTSOLUTION file
        if out_cmt_file:
            M = gcmt['moment_tensor']
            with open(out_cmt_file, 'w') as fp:
                # header line: 
                #PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
                # which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region
                fp.write(new_centroid_time.strftime(
                    'RELOC %Y %m %d %H %M %S.%f ') + \
                    '%.4f %.4f %.1f 0.0 0.0 END\n' % (evla1,evlo1,evdp1) )
                fp.write('event name:      %s\n'     % (event_id))
                fp.write('time shift:      0.0\n'                ) 
                fp.write('half duration:   %.1f\n'   % (gcmt['half_duration']))
                #fp.write('half duration:   0.0\n'    % (gcmt['half_duration']))
                fp.write('latitude:        %.4f\n'   % (evla1)   )
                fp.write('longitude:       %.4f\n'   % (evlo1)   )
                fp.write('depth:           %.4f\n'   % (evdp1)   )
                fp.write('Mrr:             %12.4e\n' % (M[0][0]) )
                fp.write('Mtt:             %12.4e\n' % (M[1][1]) )
                fp.write('Mpp:             %12.4e\n' % (M[2][2]) )
                fp.write('Mrt:             %12.4e\n' % (M[0][1]) )
                fp.write('Mrp:             %12.4e\n' % (M[0][2]) )
                fp.write('Mtp:             %12.4e\n' % (M[1][2]) )


    def plot_misfit(self, event_id, window_id, out_file=None):
        """Plot misfit for a certain event and window_id  
        """
        # CC0 map    | CC0 v.s. SNR (size ~ weight)
        #------------|-----------------
        # DTcc map   | avg. CC0            

        # check inputs
        events = self.data['events']
        if event_id not in events:
            print "[ERROR] %s does NOT exist. Exit" \
                    % (event_id)
            sys.exit()
        event = events[event_id]
        stations = event['stations']

        # get list of station,window id
        #sta_win_id_list = []
        stla_list = []
        stlo_list = []
        cc_dt_list = []
        CC0_list = []
        CCmax_list = []
        snr_list = []
        weight_list = []
        for station_id in stations:
            station = stations[station_id]
            windows = station['windows']

            # skip bad station 
            if station['stat']['code'] < 0:
                continue

            if window_id not in windows:
                continue

            window = windows[window_id]
            if window['stat']['code'] != 1:
                continue

            meta = station['meta']
            misfit = window['misfit']
            quality = window['quality']

            #sta_win_id = (station_id, window_id)
            #sta_win_id_list.append(sta_win_id)
            stla_list.append(meta['latitude'])
            stlo_list.append(meta['longitude'])
            cc_dt_list.append(misfit['CC_time_shift'])
            CC0_list.append(misfit['CC0'])
            CCmax_list.append(misfit['CCmax'])
            snr_list.append(quality['SNR'])
            weight_list.append(window['weight'])

        # get event data
        gcmt = event['gcmt']
        evla = gcmt['latitude']
        evlo = gcmt['longitude']
        M = gcmt['moment_tensor']
        Mrr = M[0][0]
        Mtt = M[1][1]
        Mpp = M[2][2]
        Mrt = M[0][1]
        Mrp = M[0][2]
        Mtp = M[1][2]
        focmec = [ Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]

        # map range
        min_lat = min(min(stla_list), evla)
        max_lat = max(max(stla_list), evla)
        lat_range = max_lat - min_lat
        min_lat -= 0.1*lat_range
        max_lat += 0.1*lat_range
        min_lon = min(min(stlo_list), evlo)
        max_lon = max(max(stlo_list), evlo)
        lon_range = max_lon - min_lon
        min_lon -= 0.1*lon_range
        max_lon += 0.1*lon_range
        #lat_true_scale = np.mean(stla_list)
        lat_0 = np.mean(stla_list)
        lon_0 = np.mean(stlo_list)
        # 
        parallels = np.arange(0.,81,10.)
        meridians = np.arange(0.,351,10.)

        # figure size
        fig = plt.figure(figsize=(11, 8.5))
        str_title = '%s %s' % (event_id, window_id)
        fig.text(0.5, 0.95, str_title, size='x-large', 
                horizontalalignment='center')

        #------ color map CC_time_shift, symbol size ~ SNR 
        ax = fig.add_axes([0.05, 0.5, 0.4, 0.35])
        ax.set_title("DT_cc (symbol_size ~ SNR)")

        m = Basemap(projection='merc', resolution='l',
                llcrnrlat=min_lat, llcrnrlon=min_lon, 
                urcrnrlat=max_lat, urcrnrlon=max_lon,
                lat_0=lat_0, lon_0=lon_0 )
        m.drawcoastlines(linewidth=0.1)
        m.drawcountries(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1])
        m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1])
        
        # CC_time_shift, SNR
        sx, sy = m(stlo_list, stla_list)
        size_list = [ 0.1 if x<0.1 else x for x in snr_list ]
        im = m.scatter(sx, sy, s=size_list, marker='o',
                c=cc_dt_list, cmap='seismic', 
                edgecolor='grey', linewidth=0.05)
        mean_amp = np.mean(cc_dt_list)
        std_amp = np.std(cc_dt_list)
        #plot_amp = abs(mean_amp)+std_amp
        plot_amp = 5.0 
        im.set_clim(-plot_amp, plot_amp)
        
        # focal mechanism
        sx, sy = m(evlo, evla)
        b = Beach(focmec, xy=(sx, sy), width=200000, linewidth=0.2, 
                facecolor='k')
        ax.add_collection(b)
        
        # colorbar
        cbar_ax = fig.add_axes([0.46, 0.575, 0.005, 0.2])
        fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=9)
        cbar_ax.set_xlabel('DT_cc(s)', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')
       
        #------ color map CC0, symbol size ~ SNR
        ax = fig.add_axes([0.05, 0.05, 0.4, 0.35])
        ax.set_title("CC0 (symbol_size ~ SNR)")

        m = Basemap(projection='merc', resolution='l',
                llcrnrlat=min_lat, llcrnrlon=min_lon, 
                urcrnrlat=max_lat, urcrnrlon=max_lon,
                lat_0=lat_0, lon_0=lon_0 )
        m.drawcoastlines(linewidth=0.1)
        m.drawcountries(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1])
        m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1])
        
        # CC0, SNR 
        sx, sy = m(stlo_list, stla_list)
        #size_list = [ 20**x for x in CCmax_list ]
        size_list = [ 0.1 if x<0.1 else x for x in snr_list ]
        im = m.scatter(sx, sy, s=size_list, marker='o',
                c=CC0_list, cmap='jet', 
                edgecolor='grey', linewidth=0.05)
        im.set_clim(0.5, 1.0)
        
        # focal mechanism
        sx, sy = m(evlo, evla)
        b = Beach(focmec, xy=(sx, sy), width=200000, linewidth=0.2, 
                facecolor='k')
        ax.add_collection(b)
 
        #add colorbar
        cbar_ax = fig.add_axes([0.46, 0.125, 0.005, 0.2])
        fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=9) 
        cbar_ax.set_xlabel('CC0', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')

        #------ SNR v.s. CC0, colored by cc_dt, size ~ weight
        ax = fig.add_axes([0.58, 0.65, 0.35, 0.2])
        im = ax.scatter(snr_list, CC0_list, marker='o', 
                s=10.*np.array(weight_list), 
                c=cc_dt_list, cmap='seismic',
                edgecolor='grey', linewidth=0.05)
        mean_amp = np.mean(cc_dt_list)
        std_amp = np.std(cc_dt_list)
        #plot_amp = abs(mean_amp)+std_amp
        plot_amp = 5.0
        im.set_clim(-plot_amp, plot_amp)
        #ax.set_xlim([min(snr_list), max(snr_list)])
        #ax.set_ylim([min(CCmax_list), max(CCmax_list)])
        ax.set_xlim([0, max(snr_list)])
        ax.set_ylim([0.3, 1.0])
        ax.set_xlabel("SNR")
        ax.set_ylabel("CC0")
        #add colorbar
        cbar_ax = fig.add_axes([0.95, 0.65, 0.005, 0.2])
        fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=9)
        cbar_ax.set_xlabel('DT_cc(s)', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')

        ##------ CC0 v.s. CCmax, colored by cc_dt
        #ax = fig.add_axes([0.58, 0.375, 0.35, 0.2])
        #im = ax.scatter(CC0_list, CCmax_list, marker='o', s=10,
        #        c=cc_dt_list, cmap='seismic',
        #        edgecolor='grey', linewidth=0.05)
        #mean_amp = np.mean(cc_dt_list)
        #std_amp = np.std(cc_dt_list)
        #plot_amp = abs(mean_amp)+std_amp
        #im.set_clim(-plot_amp, plot_amp)
        #ax.set_xlim([min(CC0_list), max(CC0_list)])
        #ax.set_ylim([min(CCmax_list), max(CCmax_list)])
        #ax.set_xlabel("CC0")
        #ax.set_ylabel("CCmax")
        ##add colorbar
        #cbar_ax = fig.add_axes([0.95, 0.375, 0.005, 0.2])
        #fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        #cbar_ax.tick_params(labelsize=9)
        #cbar_ax.set_xlabel('cc_dt(s)', fontsize=9)
        #cbar_ax.xaxis.set_label_position('top')

        ##------ cc_dt v.s. CCmax, colored by SNR
        #ax = fig.add_axes([0.58, 0.1, 0.35, 0.2])
        #im = ax.scatter(cc_dt_list, CCmax_list, marker='o', s=10, 
        #        c=snr_list, cmap='seismic',
        #        edgecolor='grey', linewidth=0.05)
        #im.set_clim(min(snr_list), max(snr_list))
        #ax.set_xlim([min(cc_dt_list), max(cc_dt_list)])
        #ax.set_ylim([min(CCmax_list), max(CCmax_list)])
        #ax.set_xlabel("cc_dt")
        #ax.set_ylabel("CCmax")
        ##add colorbar
        #cbar_ax = fig.add_axes([0.95, 0.1, 0.005, 0.2])
        #fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        #cbar_ax.tick_params(labelsize=9)
        #cbar_ax.set_xlabel('SNR(dB)', fontsize=9)
        #cbar_ax.xaxis.set_label_position('top')

        ##------ histogram of dt_cc and dt_res
        #ax1 = fig.add_axes([0.5,0.28,0.4,0.15])
        #n, bins, patches = ax1.hist(dt_cc, 50, facecolor='green', alpha=0.75)
        #amp = max(abs(dt_cc))
        #ax1.set_xlim([-amp, amp])
        #ax1.set_title('dt_cc: mean=%.2f std=%.2f' % (np.mean(dt_cc), np.std(dt_cc)))
        #ax1.tick_params(labelsize=10) 
        #
        #ax2 = fig.add_axes([0.5,0.07,0.4,0.15])
        #n, bins, patches = ax2.hist(dt_res, 50, facecolor='green', alpha=0.75)
        #amp = max(abs(dt_cc))
        #ax2.set_xlim([-amp, amp])
        #ax2.set_title('dt_res: mean=%.2f std=%.2f' % (np.mean(dt_res), np.std(dt_res)))
        #ax2.set_xlabel('dt (sec)')
        #ax2.tick_params(labelsize=10)

        #------ save figure
        if not out_file:
            out_file = '%s_%s.pdf' % (event_id, window_id)
        fig.savefig(out_file, format='pdf')
        #fig.savefig("misfit.pdf", bbox_inches='tight', format='pdf')


    def plot_seismograms(self, event_id,
            azbin=10, win=[0,100], rayp=10,
            obs_dir='obs', syn_dir='syn', syn_band_code='MX',
            syn_suffix='.sem.sac', savefig=False, out_dir='plot',
            use_STF=False, hdur=None,
            use_window=False, window_id='F.p,P',
            min_SNR=None, min_CC0=None, min_CCmax=None,
            dist_range=None):
        """ Plot seismograms for one event
            azbin:
                azimuthal bin size
        """
        event = self.data['events'][event_id]

        #====== get event data
        gcmt = event['gcmt']
        centroid_time = UTCDateTime(gcmt['centroid_time'])
        evla = gcmt['latitude']
        evlo = gcmt['longitude']
        M = gcmt['moment_tensor']
        Mrr = M[0][0]
        Mtt = M[1][1]
        Mpp = M[2][2]
        Mrt = M[0][1]
        Mrp = M[0][2]
        Mtp = M[1][2]
        focmec = [ Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]

        if use_STF:
            if not hdur:
                half_duration = gcmt['half_duration']
            else:
                half_duration = hdur
            print "half_duration = ", half_duration

        #====== get list of station,window id
        stations = event['stations']
        stla_all = []
        stlo_all = []
        for station_id in stations:
            station = stations[station_id]
            windows = station['windows']
            meta = station['meta']
            # select data 
            if station['stat']['code'] < 0:
                continue
            if use_window and (window_id not in windows):
                continue
            #
            stla_all.append(meta['latitude'])
            stlo_all.append(meta['longitude'])

        #====== calculate traveltime curve
        model = TauPyModel(model="ak135")
        dist_ttcurve = np.arange(0.0,40,0.1)
        ttcurve_p = []
        ttcurve_P = []
        ttcurve_s = []
        ttcurve_S = []
        for dist in dist_ttcurve:
            arrivals = model.get_travel_times(source_depth_in_km=gcmt['depth'],
                distance_in_degree=dist, phase_list=['p','P','s','S'])
            for arr in arrivals:
                if arr.name == 'p':
                    ttcurve_p.append((arr.distance, arr.time, arr.ray_param))
                elif arr.name == 'P':
                    ttcurve_P.append((arr.distance, arr.time, arr.ray_param))
                elif arr.name == 's':
                    ttcurve_s.append((arr.distance, arr.time, arr.ray_param))
                elif arr.name == 'S':
                    ttcurve_S.append((arr.distance, arr.time, arr.ray_param))
        
        # sort phases
        ttcurve_p = sorted(ttcurve_p, key=lambda x: x[2])
        ttcurve_P = sorted(ttcurve_P, key=lambda x: x[2])
        ttcurve_s = sorted(ttcurve_s, key=lambda x: x[2])
        ttcurve_S = sorted(ttcurve_S, key=lambda x: x[2])
        
        #====== plot map/seismograms in azimuthal bins

        #------ map configuration 
        min_lat = min(min(stla_all), evla)
        max_lat = max(max(stla_all), evla)
        lat_range = max_lat - min_lat
        min_lat -= 0.1*lat_range
        max_lat += 0.1*lat_range
        min_lon = min(min(stlo_all), evlo)
        max_lon = max(max(stlo_all), evlo)
        lon_range = max_lon - min_lon
        min_lon -= 0.1*lon_range
        max_lon += 0.1*lon_range
        #lat_true_scale = np.mean(stla_list)
        lat_0 = np.mean(stla_all)
        lon_0 = np.mean(stlo_all)
        # 
        parallels = np.arange(0.,81,10.)
        meridians = np.arange(0.,351,10.)

        #------ plot each azimuthal bin
        for az in np.arange(0, 360, azbin):
            azmin = az
            azmax = az + azbin
            # read available stations within azbin
            data_azbin = {}
            for station_id in stations:
                station = stations[station_id]
                meta = station['meta']
                azimuth = meta['azimuth']
                windows = station['windows']
                dist = meta['dist_degree']
                if station['stat']['code'] < 0:
                    continue
                if dist_range:
                    if dist < min(dist_range) or dist > max(dist_range):
                        continue
                if azimuth<azmin or azimuth>azmax:
                    continue
                if use_window: 
                    if (window_id not in windows):
                        continue
                    window = windows[window_id]
                    quality = window['quality']
                    misfit = window['misfit']
                    if window['stat']['code'] != 1:
                        continue
                    if min_SNR and quality['SNR']<min_SNR:
                        continue
                    if min_CC0 and misfit['CC0']<min_CC0:
                        continue
                    if min_CCmax and misfit['CCmax']<min_CCmax:
                        continue
                try:
                    syn_st, obs_st = self.read_seismograms(
                        event_id, station_id,
                        obs_dir=obs_dir, syn_dir=syn_dir,
                        syn_band_code=syn_band_code, syn_suffix=syn_suffix)
                except Exception as e:
                    print str(e)
                    continue
                # syn/obs data arrays
                syn_npts = syn_st[0].stats.npts
                syn_delta = syn_st[0].stats.delta
                syn_nyq = 0.5 / syn_delta
                syn_ENZ = np.zeros((3, syn_npts))
                obs_ENZ = np.zeros((3, syn_npts))
                for i in range(3):
                    syn_ENZ[i,:] = syn_st[i].data
                    obs_ENZ[i,:] = obs_st[i].data

                # desgin filter
                filter_param = station['filter']
                filter_type = filter_param['type']
                filter_order = filter_param['order']
                filter_freqlim = filter_param['freqlim']
                filter_b, filter_a = signal.butter(filter_order,
                        np.array(filter_freqlim)/syn_nyq, btype='band')
                # filter obs: F * d
                obs_ENZ[:,:] = signal.filtfilt(filter_b, filter_a, obs_ENZ)
                # filter syn: F * S * u
                syn_ENZ[:,:] = signal.filtfilt(filter_b, filter_a, syn_ENZ)
                if use_STF:
                    f = np.fft.rfftfreq(syn_npts, d=syn_delta)
                    F_src = np.sinc(f * half_duration)**2
                    syn_ENZ[:,:] = np.fft.irfft(F_src*np.fft.rfft(syn_ENZ), syn_npts)

                # rotate EN -> TR (TRZ: right-hand convention)
                Raz = (meta['back_azimuth'] + 180.0) % 360.0
                sin_Raz = np.sin(np.deg2rad(Raz))
                cos_Raz = np.cos(np.deg2rad(Raz))
                proj_matrix = [ [sin_Raz,  cos_Raz],
                                [cos_Raz, -sin_Raz] ]
                syn_ENZ[0:2,:] = np.dot(proj_matrix, syn_ENZ[0:2,:])
                obs_ENZ[0:2,:] = np.dot(proj_matrix, obs_ENZ[0:2,:])
                syn_RTZ = syn_ENZ
                obs_RTZ = obs_ENZ
                # record results
                starttime = syn_st[0].stats.starttime
                dt = syn_st[0].stats.delta
                npts = syn_st[0].stats.npts
                times = dt * np.arange(npts)
                data_azbin[station_id] = {'meta':meta, 'starttime':starttime,
                        'times':times, 'syn':syn_RTZ, 'obs':obs_RTZ}

            if not data_azbin: continue

            #====== create figure and axes
            fig = plt.figure(figsize=(8.5, 11)) # US letter
            str_title = '{:s} (win: {:s}, az: {:04.1f}~{:04.1f})'.format(
                    event_id, window_id, azmin, azmax)
            fig.text(0.5, 0.95, str_title, size='x-large', horizontalalignment='center')

            #------ station map
            ax_origin = [0.3, 0.74]
            ax_size = [0.4, 0.2]
            ax_map = fig.add_axes(ax_origin + ax_size)
            m = Basemap(projection='merc', resolution='l',
                    llcrnrlat=min_lat, llcrnrlon=min_lon, 
                    urcrnrlat=max_lat, urcrnrlon=max_lon,
                    lat_0=lat_0, lon_0=lon_0 )
            m.drawcoastlines(linewidth=0.1)
            m.drawcountries(linewidth=0.1)
            m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1], 
                    fontsize=10, fmt='%3.0f')
            m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1], 
                    fontsize=10, fmt='%3.0f')
            sx, sy = m(stlo_all, stla_all)
            m.scatter(sx, sy, s=10, marker='^', facecolor='blue', edgecolor='')
            # plot stations inside the bin
            stla = [data_azbin[x]['meta']['latitude'] for x in data_azbin]
            stlo = [data_azbin[x]['meta']['longitude'] for x in data_azbin]
            sx, sy = m(stlo, stla)
            m.scatter(sx, sy, s=10, marker='^', facecolor='red', edgecolor='')

            # focal mechanism
            sx, sy = m(evlo, evla)
            b = Beach(focmec, xy=(sx, sy), width=400000, linewidth=0.2, facecolor='r')
            ax_map.add_collection(b)
 
            #------ plot waveforms 
            ax_RTZ = []
            for i in range(3):
                ax_origin = [0.07+0.3*i, 0.05]
                ax_size = [0.25, 0.65]
                ax_RTZ.append(fig.add_axes(ax_origin + ax_size))

            y = [ x['meta']['dist_degree'] for x in data_azbin.itervalues() ]
            ny = len(y)
            dy = 0.5*(max(y)-min(y)+1)/ny
            if dist_range:
                plot_ymax = max(dist_range) + 2*dy
                plot_ymin = min(dist_range) - 2*dy
            else:
                plot_ymax = max(y) + 2*dy
                plot_ymin = min(y) - 2*dy
        
            #plot traveltime curves 
            for i in range(3):
                ax = ax_RTZ[i]
                ax.plot([x[1]-rayp*x[0] for x in ttcurve_p], [x[0] for x in ttcurve_p], 'b-', linewidth=0.2)
                ax.plot([x[1]-rayp*x[0] for x in ttcurve_P], [x[0] for x in ttcurve_P], 'b-', linewidth=0.2)
                ax.plot([x[1]-rayp*x[0] for x in ttcurve_s], [x[0] for x in ttcurve_s], 'c-', linewidth=0.2)
                ax.plot([x[1]-rayp*x[0] for x in ttcurve_S], [x[0] for x in ttcurve_S], 'c-', linewidth=0.2)
 
            cmp_names = ['R', 'T', 'Z']
            for station_id in data_azbin:
                sta = data_azbin[station_id]
                meta = sta['meta']
                dist_degree = meta['dist_degree']
                reduced_time = dist_degree * rayp
                # time of first sample referred to centroid time 
                t0 = sta['starttime'] - centroid_time
                # time of samples referred to centroid time
                t = sta['times'] + t0
                plot_t0 = win[0] + reduced_time
                plot_t1 = win[1] + reduced_time
                idx = (t > plot_t0) & (t < plot_t1)
                    
                t_plot = t[idx] - reduced_time
                obs_RTZ = sta['obs']
                syn_RTZ = sta['syn']

                windows = stations[station_id]['windows']
                if use_window and (window_id in windows):
                    window = windows[window_id]
                    quality = window['quality']
                    A_obs = quality['A_obs']
                    A_syn = quality['A_syn']
                    win_starttime = UTCDateTime(window['starttime'])
                    win_endtime = UTCDateTime(window['endtime'])
                    win_t0 = win_starttime - centroid_time - reduced_time
                    win_t1 = win_endtime - centroid_time - reduced_time
                else:
                    A_obs = np.sqrt(np.max(np.sum(obs_RTZ[:,idx]**2, axis=0)))
                    A_syn = np.sqrt(np.max(np.sum(syn_RTZ[:,idx]**2, axis=0)))

                for i in range(3):
                    #normalize data
                    obs = obs_RTZ[i, idx]
                    obs = dy*obs/A_obs
                    syn = syn_RTZ[i, idx]
                    syn = dy*syn/A_syn

                    ax = ax_RTZ[i]
                    ax.plot(t_plot, obs+dist_degree, 'k-', linewidth=0.5)
                    ax.plot(t_plot, syn+dist_degree, 'r-', linewidth=0.5)

                    # annotatate time window
                    if use_window:
                        ax.plot(win_t0, dist_degree, 'k|', markersize=8)
                        ax.plot(win_t1, dist_degree, 'k|', markersize=8)
                        misfit = window['misfit']
                        # CC0
                        if i == 0:
                            ax.text(win[1], dist_degree, ' %.3f' % (misfit['CC0']), 
                                    verticalalignment='center', fontsize=7)
                        # window weight
                        if i == 1:
                            ax.text(win[1], dist_degree, ' %.1f' % (window['weight']), 
                                    verticalalignment='center', fontsize=7)

                    #annotate station names 
                    if i == 2:
                        #str_annot = '%.3f,%.1f,%s' % (
                        #        misfit['CC0'], window['weight'], station_id)
                        ax.text(win[1], dist_degree, ' '+station_id, \
                                verticalalignment='center', fontsize=7)

                #for i in range(3):
            #for sta_id in data:
        
            # control axes limits and lables, annotation
            for i in range(3):
                # axis 
                ax.set_xlim(win[0], win[1])
                ax.set_ylim(plot_ymin, plot_ymax)
                ax.set_title(cmp_names[i])
                ax.set_xlabel('t - {:.1f}*dist (s)'.format(rayp))
                ax.tick_params(axis='both',labelsize=10)
                # ylabel 
                if i == 0:
                    ax.set_ylabel('dist (deg)')
                else:
                    ax.set_yticklabels([])

            # save figures
            if savefig:
                if use_window:
                    out_file = '%s/%s_az_%03d_%03d_%s.pdf' % (
                            out_dir, event_id, azmin, azmax, window_id)
                else:
                    out_file = '%s/%s_az_%03d_%03d.pdf' % (
                            out_dir, event_id, azmin, azmax)
                plt.savefig(out_file, format='pdf')
            else:
                plt.show()

            plt.close(fig)
        
        # END