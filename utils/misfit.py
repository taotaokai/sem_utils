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


#====== utility functions
def is_equal(lst):
    return not lst or [lst[0]]*len(lst) == lst


def cosine_taper(x, xc):
    """cosine taper at two ends stop,pass, pass,stop
        xc: (array-like)
            stop,pass[,pass,stop]
        x: array-like
            sample points
    """
    nc = len(xc)
    x = np.array(x)
    y = np.ones(len(x))

    if nc == 2: # sided taper
        l = xc[1] - xc[0] # taper width
        if l == 0:
            raise ValueError('pass and stop values cannot be the same.')
        elif l > 0: # tapered at left side
            idx = x < xc[0]; y[idx] = 0.0
            idx = (xc[0] <= x) & (x <= xc[1])
            y[idx] = 0.5 - 0.5*np.cos(np.pi*(x[idx] - xc[0])/l)
            idx = x > xc[1]; y[idx] = 1.0
        else: # tapered at right side
            idx = x > xc[0]; y[idx] = 0.0
            idx = (xc[1] <= x) & (x <= xc[0])
            y[idx] = 0.5 + 0.5*np.cos(np.pi*(x[idx] - xc[1])/l)
            idx = x < xc[1]; y[idx] = 1.0
    elif nc == 4: # two sided taper
        if not (xc[0]<xc[1]<xc[2]<xc[3]):
            raise ValueError('4 cutoff values must be in increasing order.')
        else:
            idx = x <= xc[0]; y[idx] = 0.0

            idx = (xc[0] < x) & (x < xc[1])
            y[idx] = 0.5 - 0.5*np.cos(np.pi*(x[idx] - xc[0])/(xc[1]-xc[0]))
            
            idx = (xc[1] <= x) & (x <= xc[2]); y[idx] = 1.0

            idx = (xc[2] < x) & (x < xc[3])
            y[idx] = 0.5 + 0.5*np.cos(np.pi*(x[idx] - xc[2])/(xc[3]-xc[2]))

            idx = x > xc[3]; y[idx] = 0.0
    else:
        raise ValueError('number of cutoff values must be either 2 or 4.')

    return y


def taper_window(npts, ttype="cosine", ratio=0.1):
    """taper window of npts points
        ratio: decimal percentage at one end (default: 0.1)
    """
    if not 0.0 < ratio <= 0.5:
        raise ValueError("taper ratio must be between 0.0 and 0.5")

    xc = [0.0, npts*ratio, npts*(1.0-ratio), npts]
    x = np.arange(npts)
    if ttype == "cosine":
        taper = cosine_taper(x, xc)
    else:
        raise Exception("unrecognized taper type" % (ttype))

    return taper

#======
class Misfit(object):
    """Class managing all misfit windows

    data structure:

    Misfit: {
    |   'iteration': {iter, ...}, 
    |   'events': {
    |   *   <event_id>: {
    |   |   |   'gcmt': {code, lat,lon,depth, ...}, 
    |   |   |   'stations': {
    |   |   |   *   <station_id(net.sta.loc)>: {
    |   |   |   |   |   'stat': {code, msg}, # code<0: problematic, code=0: OK
    |   |   |   |   |   'meta': {lat,lon,ele,...,channels:[{code,az,dip,...},...]},
    |   |   |   |   |   'filter': {type, freqlim(list)},
    |   |   |   |   |   'taper': {type, ratio},
    |   |   |   |   |   'noise_window': {starttime, endtime}, # pre-event time window
    |   |   |   |   |   'windows': {
    |   |   |   |   |   *   <window_id(cmp.pha)>: {
    |   |   |   |   |   |   |   component, azimuth, dip, starttime, endtime, weight,
    |   |   |   |   |   |   |   'stat': {code, msg}, #code<0: problematic, code=0: not measured, code=1: measured
    |   |   |   |   |   |   |   'phase': {name, ttime, takeoff_angle, ray_param},
    |   |   |   |   |   |   |   'quality': {A_obs, A_noise, SNR},
    |   |   |   |   |   |   |   'misfit': {cc_0, cc_max, cc_time_shift,
    |   |   |   |   |   |   |              ar_0, ar_max } },
    |   |   |   |   |   *   <window_id>: { },
    |   |   |   |   |   |   ...
    |   |   |   *   <station_id>: { ... },
    |   |   |   |   ...
    |   *   <event_id>: { ... },
    |   |   ...
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
            noise_window=('ttp',[-150,-40]),
            filter_param=('cosine', [0.009, 0.011, 0.06, 0.10]),
            taper_param=('cosine', 0.1),
            event_id_list=None, station_id_list=None, update=False):
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

        # noise window param
        noise_phases = noise_window[0].split(',')
        noise_begin = noise_window[1][0]
        noise_end = noise_window[1][1]

        # filter/taper parameters
        filter_param = {'type': filter_param[0], 'freqlim': filter_param[1]}
        taper_param = {'type': taper_param[0], 'ratio': taper_param[1]}

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
                station['filter'] = filter_param
                station['taper'] = taper_param

                # setup noise window
                # make noise window
                arrivals = taup_model.get_travel_times(
                        source_depth_in_km=gcmt['depth'],
                        distance_in_degree=dist_degree,
                        phase_list=noise_phases)
                if arrivals:
                    arr = arrivals[0]
                    noise_starttime = str(centroid_time + arr.time + noise_begin)
                    noise_endtime = str(centroid_time + arr.time + noise_end)
                else:
                    print "[WARNNING] %s:%s:%s phase(s) not found, skip." \
                            % (event_id, station_id, noise_window[0])
                    continue
                noise_window = {'starttime':noise_starttime, 'endtime':noise_endtime}
                station['noise_window'] = noise_window

                # setup singal windows
                windows = station['windows']
                for win in window_list:
                    comp = win[0]
                    phase = win[1]
                    signal_begin = float(win[2][0])
                    signal_end = float(win[2][1])
                    window_id = "%s.%s" % (comp, phase)
                    if (window_id in windows) and (not update):
                        print "[WARNING] %s:%s:%s already exists, skip" \
                                % (event_id, station_id, window_id)
                        continue
                    if comp == 'Z': # vertcal component
                        cmpaz = 0.0 
                        cmpdip = -90.0
                    elif comp == 'R': # radial component
                        cmpaz = (baz + 180.0)%360.0
                        cmpdip = 0.0
                    elif comp == 'T': # tangential component
                        cmpaz = (baz + 90.0)%360.0
                        cmpdip = 0.0
                    elif comp == 'H': # horizontal vector
                        cmpaz = float('nan')
                        cmpdip = float('nan')
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

                    # add/update window
                    # only use first arrival of smallest traveltime
                    arr = arrivals[0]
                    signal_starttime = str(centroid_time + arr.time + signal_begin)
                    signal_endtime = str(centroid_time + arr.time + signal_end)
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
                    if window_id not in windows:
                        windows[window_id] = window
                    elif update:
                        windows[window_id].update(window)
                        windows[window_id]['stat']['msg'] = "updated"

                #end for winpar in winpar_list:
            #end for station_id in station_id_list:
        #end for event_id in event_id_list:


    def read_seismograms_for_one_station(self, event_id, station_id, 
            obs_dir='obs', syn_dir='syn', syn_band_code='MX', 
            syn_suffix='.sem.sac', filter_padlen=50,
            use_prefilt=False, prefilt=('cosine', [0.13, 0.07]),
            use_STF=False, STF=('triangle', 0.0)):
        """ Get 3-component seismograms recorded at one station for one event.
            filter_padlen:
                 length of zero-padding at the end of syn data before filtering.
            use_prefilt, prefilt:
                lowpass filter which is set as the source time function in 
                forward simulations to suppress high-frequency numerical noise.
                if selected, this is applied to obs data.
            use_STF, STF:
                if selected, source time function is applied to syn data.
        """
        syn_orientation_codes = ['E', 'N', 'Z']
        events = self.data['events']
        event = events[event_id]
        stations = event['stations']
        station = stations[station_id]
 
        #------ get file paths of obs, syn seismograms
        meta = station['meta']
        channels = meta['channels']
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
        syn_delta = tr.stats.delta
        syn_sampling_rate = tr.stats.sampling_rate
        syn_npts = tr.stats.npts
        syn_times = syn_delta*np.arange(syn_npts)

        #------ filter parameters
        filter_param = station['filter']
        filter_type = filter_param['type']
        filter_freqlim = filter_param['freqlim']

        #------ filter syn seismograms
        npad = int(filter_padlen/syn_delta)
        f = np.fft.rfftfreq(syn_npts + npad, d=syn_delta)
        H = cosine_taper(f, filter_freqlim)
        if use_STF:
            # source spectrum: isosceles triangle
            half_duration = STF[1]
            src_spectrum = np.sinc(f * half_duration)**2
            H *= src_spectrum
        syn_ENZ = np.zeros((3, syn_npts))
        for i in range(3):
            x = np.fft.irfft(H * np.fft.rfft(np.concatenate(
                (syn_st[i].data, np.zeros(npad)))), syn_npts+npad)
            syn_ENZ[i,:] = x[0:syn_npts]

        #------ process obs seismograms
        # noise window
        noise_starttime = UTCDateTime(station['noise_window']['starttime'])
        noise_endtime = UTCDateTime(station['noise_window']['endtime'])
        noise_window_len = noise_endtime - noise_starttime
        noise_npts = int(noise_window_len / syn_delta)
        noise_times = syn_delta * np.arange(noise_npts)
        # filter obs, then cut into signal/noise windows
        obs_ENZ = np.zeros((3, syn_npts))
        noise_ENZ = np.zeros((3, noise_npts))
        for i in range(3):
            tr = obs_st[i]
            # fitler
            x = signal.detrend(tr.data, type='linear')
            obs_npts = tr.stats.npts
            obs_delta = tr.stats.delta
            #npad = int(filter_padlen/obs_delta)
            f = np.fft.rfftfreq(obs_npts, d=obs_delta)
            H = cosine_taper(f, filter_freqlim)
            if use_prefilt:
                H *= cosine_taper(f, prefilt[1])
            x = np.fft.irfft(H * np.fft.rfft(x), obs_npts)
            # interpolate obs into the same time samples of syn
            obs_ENZ[i,:] = lanczos_interp1(x, tr.stats.delta,
                    syn_times+(syn_starttime-tr.stats.starttime), na=20)
            # interpolate obs into noise window
            noise_ENZ[i,:] = lanczos_interp1(x, tr.stats.delta,
                    noise_times+(noise_starttime-tr.stats.starttime), na=20)
        # roate obs to ENZ
        #   projection matrix: obs = proj * ENZ => ZNE = inv(proj) * obs
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
        noise_ENZ = np.dot(inv_proj, noise_ENZ)

        return syn_starttime, syn_times, syn_delta, syn_ENZ, obs_ENZ, \
                noise_starttime, noise_ENZ


    def measure_windows_for_one_event(self, event_id,
            obs_dir='obs', syn_dir='syn', syn_band_code='MX', 
            syn_suffix='.sem.sac',
            use_prefilt=False, prefilt=(0.13, 0.07),
            use_STF=False,     STF=('triangle', 0.0),
            cc_delta=0.01, output_adj=False, adj_dir='adj', 
            adj_window_id_list=['F.p,P','F.s,S'],
            weight_param={'SNR':[5,10], 'cc_max':[0.6,0.8], 'cc_0':[0.5,0.7]},
            update=False):
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
            print "[ERROR] %s does NOT exist."  % (event_id)
            sys.exit()
        event = events[event_id]
        
        # process each station
        stations = event['stations']
        for station_id in stations:
            station = stations[station_id]
            # window taper parameters
            taper_param = station['taper']
            taper_type = taper_param['type']
            taper_ratio = taper_param['ratio']

            # read seismograms
            syn_startime, syn_times, syn_detla, syn_ENZ, obs_ENZ, \
                    noise_starttime, noise_ENZ = \
                self.read_seismograms_for_one_station(event_id,
                        station_id, obs_dir=obs_dir, syn_dir=syn_dir,
                        syn_band_code=syn_band_code,
                        syn_suffix=syn_suffix,
                        filter_padlen=filter_padlen, 
                        use_prefilt=use_prefilt, prefilt=prefilt, 
                        use_STF=use_STF, STF=STF)
            syn_sampling_rate = 1.0/syn_delta

            # process each signal window
            windows = station['windows']
            if output_adj:
                adj_ENZ = np.zeros((3, syn_npts))
            taper = np.zeros(syn_npts)

            for window_id in windows:
                window = windows[window_id]

                comp = window['component']
                cmpaz = window['azimuth']
                cmpdip = window['dip']
                window_starttime = UTCDateTime(window['starttime'])
                window_endtime = UTCDateTime(window['endtime'])
                window_len = window_endtime - window_starttime

                # make taper window
                win_b = window_starttime - syn_starttime
                win_e = window_endtime - syn_starttime
                win_ib = int(win_b * syn_sampling_rate)
                win_ie = int(win_e * syn_sampling_rate) + 1
                if win_ib < 0: win_ib = 0
                if win_ie > syn_npts: win_ie = syn_npts
                taper[:] = 0.0
                taper[win_ib:win_ie] = taper_window(win_ie-win_ib, 
                        taper_type, taper_ratio)
 
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
                noise_ENZ_win = np.dot(proj_matrix, noise_ENZ)
                obs_ENZ_win = np.dot(proj_matrix, obs_ENZ) * taper
                syn_ENZ_win = np.dot(proj_matrix, syn_ENZ) * taper

                # measure SNR
                A_syn = np.sqrt(np.max(np.sum(syn_ENZ_win**2, axis=0)))
                A_obs = np.sqrt(np.max(np.sum(obs_ENZ_win**2, axis=0)))
                A_noise =  np.sqrt(np.max(np.sum(noise_ENZ_win**2, axis=0)))
                if A_obs==0: # data file may not containe the time window
                    print '[WARNING] %s:%s:%s empty obs trace, SKIP.' \
                            % (event_id, station_id, window_id)
                    continue
                if A_noise==0: # could occure when the data begin time is too close to the first arrival
                    print '[WARNING] %s:%s:%s empty noise trace, SKIP.' \
                            % (event_id, station_id, window_id)
                    continue
                snr = 20.0 * np.log10(A_obs/A_noise)
 
                # measure misfit between observed and synthetic seismograms
                # cross-correlation
                obs_norm2 = np.sum(obs_ENZ_win**2)
                obs_norm = np.sqrt(obs_norm2)
                syn_norm = np.sqrt(np.sum(syn_ENZ_win**2))
                cc = np.zeros(2*syn_npts-1)
                for i in range(3):
                    # NOTE the order (obs,syn) is important. The positive time on 
                    #   CC means shifting syn in the positive time direction
                    cc += signal.fftconvolve(
                            obs_ENZ_win[i,:], syn_ENZ_win[i,::-1], 'full')
                cc /= obs_norm * syn_norm
                # zero-lag normalized correlation
                cc_0 = cc[syn_npts]
                ar_0 = cc_0 * syn_norm / obs_norm # amplitude ratio syn/obs 
                # tshift>0: synthetic is shifted along the positive time direction
                cc_max_time_shift = window_len/2.0 #TODO: more reasonable choice?
                ncc = int(cc_max_time_shift / cc_delta)
                cc_times = np.arange(-ncc,ncc+1) * cc_delta
                # interpolate cc to finer time samples
                if syn_delta < cc_delta:
                    print '[WARNING] syn_delta(%f) < cc_time_step(%f)' \
                            % (syn_delta, cc_delta)
                ti = (syn_npts-1)*syn_delta + cc_times # -(npts-1)*dt: begin time in cc
                cci = lanczos_interp1(cc, syn_delta, ti, na=20)
                # time shift at the maximum correlation
                imax = np.argmax(cci)
                cc_time_shift = cc_times[imax]
                cc_max = cci[imax]
                ar_max = cc_max * syn_norm / obs_norm # amplitude ratio: syn/obs

                # window weighting based on SNR and misfit
                weight = cosine_taper(weight_param['SNR'], snr) \
                        * cosine_taper(weight_param['cc_max'], cc_max) \
                        * cosine_taper(weight_param['cc_0'], cc_0)

                # make adjoint source
                # adjoint source type = normalized cc
                # adj = weight * (obs - A0*syn) / (obs_norm*syn_norm)
                A0 = cc_0 * obs_norm / syn_norm # amplitude raito: obs/syn
                if output_adj and (window_id in adj_window_id_list):
                    # adjoint source of the current window
                    adj_win = signal.filtfilt(filter_b, filter_a,
                            taper*(obs_ENZ_win - A0*syn_ENZ_win))/obs_norm/syn_norm
                    if use_STF:
                        adj_win = np.fft.irfft(
                            np.fft.rfft(adj_win) * np.conjugate(src_spectrum))
                    adj_ENZ += weight * adj_win

                # form measurment results
                quality_dict = {
                        'A_obs': A_obs, 'A_syn': A_syn, 'A_noise': A_noise,
                        'SNR': snr}
                misfit_dict = {
                        'cc_0': cc_0, 'cc_max': cc_max,
                        'ar_0': ar_0, 'ar_max': ar_max,
                        'cc_time_shift': cc_time_shift }

                # add/update measure results
                if 'misfit' not in window \
                        or not window['misfit']:
                    window['quality'] = quality_dict
                    window['misfit'] = misfit_dict
                    window['weight'] = weight
                    window['stat'] = {'code': 1, 'msg': "measured"}
                elif update:
                    window['quality'].update(quality_dict)
                    window['misfit'].update(misfit_dict)
                    window['weight'] = weight
                    window['stat'] = {'code': 1, 'msg': "updated"}
                else:
                    continue

                # DEBUG
                #t = syn_times
                #noise_b = noise_starttime - syn_starttime
                #for i in range(3):
                #    plt.subplot(411+i)
                #    if i == 0:
                #        plt.title('%s dt %.2f ccmax %.3f armax %.3f ' \
                #                'cc0 %.1f ar0 %.3f \nAobs %g Anoise %g snr %.1f' \
                #                % (station_id, cc_time_shift, 
                #                   cc_max, ar_max, cc_0, ar_0, A_obs, A_noise, snr) )
                #    A_plot = max([A_obs, A_syn])

                #    plt.plot(t, obs_ENZ[i,:]/A_plot, 'k', linewidth=0.2)
                #    plt.plot(noise_times+noise_b, noise_ENZ_win[i,:]/A_plot, 'b', linewidth=0.5)
                #    plt.plot(t, syn_ENZ[i,:]/A_plot, 'r', linewidth=0.2)
                #    plt.plot(t[win_ib:win_ie], obs_ENZ_win[i,win_ib:win_ie]/A_plot,
                #            'k', linewidth=1.0)
                #    plt.plot(t[win_ib:win_ie], syn_ENZ_win[i,win_ib:win_ie]/A_plot,
                #            'r', linewidth=1.0)
                #    plt.ylim((-1.5, 1.5))
                #    plt.xlim((noise_b, max(t)))
                #    plt.ylabel(syn_orientation_codes[i])

                #plt.subplot(414)
                #plt.plot(cc_times, cci, 'k-')
                #plt.xlim((min(cc_times), max(cc_times)))
                #plt.ylabel(window_id)
            
                #plt.show()

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

        #end for station_id in station_id_list:
    #end def setup_windows(self,


    def relocate_1d(self, event_id, window_id_list=['Z.p,P'],
        min_SNR=10.0, min_cc_0=0.5, min_cc_max=0.7, 
        max_cc_time_shift=8.0, fix_depth=False, out_cmt_file=None):
        """relocate event using ray path in reference earth model
        """
        # check inputs
        events = self.data['events']
        if event_id not in events:
            print "[ERROR] %s does NOT exist. Exit" \
                    % (event_id)
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
                if window['quality']['SNR'] < min_SNR or \
                        misfit['cc_0'] < min_cc_0 or \
                        misfit['cc_max'] < min_cc_max or\
                        abs(misfit['cc_time_shift']) > max_cc_time_shift:
                    continue

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
            G[i,:] = [-pn, -pe, -pd, 1.0] # -p: from receiver to source
            dt_cc[i] = misfit['cc_time_shift']

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
        reloc_dict = {
                'window_id_list': window_id_list,
                'min_SNR': min_SNR,
                'min_cc_0': min_cc_0,
                'min_cc_max': min_cc_max,
                'max_cc_time_shift': max_cc_time_shift,
                'singular_value': sigval.tolist(),
                'dm': {'dNorth':dm[0], 'dEast':dm[1], 'dDepth':dm[2], 
                    'dT':dm[3]},
                'old_gcmt': gcmt,
                'new_location': {'latitude':evla1, 'longitude':evlo1, 
                    'depth':evdp1},
                'data': {'num':n, 'mean':np.mean(dt_cc), 'std':np.std(dt_cc)},
                'residual': {'mean':np.mean(dt_res), 'std':np.std(dt_res)} }

        # make new CMTSOLUTION file
        new_centroid_time = UTCDateTime(gcmt['centroid_time']) + ev_dt
        if out_cmt_file:
            M = gcmt['moment_tensor']
            with open(out_cmt_file, 'w') as fp:
                # header line: 
                #PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
                # which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region
                fp.write(new_centroid_time.strftime(
                    'RELOC %Y %m %d %H %M %S.%f ') + \
                    '%.4f %.4f %.1f 0.0 0.0 REGION\n' % (evla1,evlo1,evdp1) )
                fp.write('event name:      %s\n'     % (event_id))
                fp.write('time shift:      0.0\n'                ) 
                fp.write('half duration:   %.1f\n'   % (gcmt['half_duration']))
                fp.write('latitude:        %.4f\n'   % (evla1)   )
                fp.write('longitude:       %.4f\n'   % (evlo1)   )
                fp.write('depth:           %.4f\n'   % (evdp1)   )
                fp.write('Mrr:             %12.4e\n' % (M[0][0]) )
                fp.write('Mtt:             %12.4e\n' % (M[1][1]) )
                fp.write('Mpp:             %12.4e\n' % (M[2][2]) )
                fp.write('Mrt:             %12.4e\n' % (M[0][1]) )
                fp.write('Mrp:             %12.4e\n' % (M[0][2]) )
                fp.write('Mtp:             %12.4e\n' % (M[1][2]) )

        return reloc_dict


    def plot_misfit(self, event_id, window_id, out_file=None):
        """Plot misfit for a certain event and window_id  
        """
        # dt_cc       | cc_0/cc_max V.S. dt_cc 
        #-------------|-----------------------
        # hist?       | cc_0/cc_max VS SNR

        # check inputs
        events = self.data['events']
        if event_id not in events:
            print "[ERROR] %s does NOT exist. Exit" \
                    % (event_id)
            sys.exit()
        event = events[event_id]
        stations = event['stations']

        # get list of data 
        #sta_win_id_list = []
        stla_list = []
        stlo_list = []
        cc_dt_list = []
        cc_0_list = []
        cc_max_list = []
        snr_list = []
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
            cc_dt_list.append(misfit['cc_time_shift'])
            cc_0_list.append(misfit['cc_0'])
            cc_max_list.append(misfit['cc_max'])
            snr_list.append(quality['SNR'])

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
        min_lat = min(min(stla_list), evla) - 1.0
        max_lat = max(max(stla_list), evla) + 1.0
        min_lon = min(min(stlo_list), evlo) - 1.0
        max_lon = max(max(stlo_list), evlo) + 1.0
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

        #------ cc_time_shift + SNR 
        ax = fig.add_axes([0.05, 0.5, 0.4, 0.35])
        ax.set_title("cc_time_shift (symbol_size ~ SNR)")

        m = Basemap(projection='merc', resolution='l',
                llcrnrlat=min_lat, llcrnrlon=min_lon, 
                urcrnrlat=max_lat, urcrnrlon=max_lon,
                lat_0=lat_0, lon_0=lon_0 )
        m.drawcoastlines(linewidth=0.1)
        m.drawcountries(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1])
        m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1])
        
        # cc_time_shift, SNR
        sx, sy = m(stlo_list, stla_list)
        size_list = [ 0.01 if x<0.01 else x for x in snr_list ]
        im = m.scatter(sx, sy, s=size_list, marker='o',
                c=cc_dt_list, cmap='seismic', 
                edgecolor='grey', linewidth=0.05)
        mean_amp = np.mean(cc_dt_list)
        std_amp = np.std(cc_dt_list)
        plot_amp = abs(mean_amp)+std_amp
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
        cbar_ax.set_xlabel(' cc_dt(s)', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')
       
        #------ color: cc_max, size: SNR
        ax = fig.add_axes([0.05, 0.05, 0.4, 0.35])
        ax.set_title("cc_max (symbol_size ~ SNR)")

        m = Basemap(projection='merc', resolution='l',
                llcrnrlat=min_lat, llcrnrlon=min_lon, 
                urcrnrlat=max_lat, urcrnrlon=max_lon,
                lat_0=lat_0, lon_0=lon_0 )
        m.drawcoastlines(linewidth=0.1)
        m.drawcountries(linewidth=0.1)
        m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1])
        m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1])
        
        # cc_max, SNR 
        sx, sy = m(stlo_list, stla_list)
        #size_list = [ 20**x for x in cc_max_list ]
        size_list = [ 0.01 if x<0.01 else x for x in snr_list ]
        im = m.scatter(sx, sy, s=size_list, marker='o',
                c=cc_max_list, cmap='seismic', 
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
        cbar_ax.set_xlabel(' cc_max', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')

        #------ SNR v.s. cc_max, colored by cc_dt
        ax = fig.add_axes([0.58, 0.65, 0.35, 0.2])
        im = ax.scatter(snr_list, cc_max_list, marker='o', s=10,
                c=cc_dt_list, cmap='seismic',
                edgecolor='grey', linewidth=0.05)
        mean_amp = np.mean(cc_dt_list)
        std_amp = np.std(cc_dt_list)
        plot_amp = abs(mean_amp)+std_amp
        im.set_clim(-plot_amp, plot_amp)
        ax.set_xlim([min(snr_list), max(snr_list)])
        ax.set_ylim([min(cc_max_list), max(cc_max_list)])
        ax.set_xlabel(" SNR")
        ax.set_ylabel("cc_max")
        #add colorbar
        cbar_ax = fig.add_axes([0.95, 0.65, 0.005, 0.2])
        fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=9)
        cbar_ax.set_xlabel('cc_dt(s)', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')

        #------ cc_0 v.s. cc_max, colored by cc_dt
        ax = fig.add_axes([0.58, 0.375, 0.35, 0.2])
        im = ax.scatter(cc_0_list, cc_max_list, marker='o', s=10,
                c=cc_dt_list, cmap='seismic',
                edgecolor='grey', linewidth=0.05)
        mean_amp = np.mean(cc_dt_list)
        std_amp = np.std(cc_dt_list)
        plot_amp = abs(mean_amp)+std_amp
        im.set_clim(-plot_amp, plot_amp)
        ax.set_xlim([min(cc_0_list), max(cc_0_list)])
        ax.set_ylim([min(cc_max_list), max(cc_max_list)])
        ax.set_xlabel("cc_0")
        ax.set_ylabel("cc_max")
        #add colorbar
        cbar_ax = fig.add_axes([0.95, 0.375, 0.005, 0.2])
        fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=9)
        cbar_ax.set_xlabel('cc_dt(s)', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')

        #------ cc_dt v.s. cc_max, colored by SNR
        ax = fig.add_axes([0.58, 0.1, 0.35, 0.2])
        im = ax.scatter(cc_dt_list, cc_max_list, marker='o', s=10, 
                c=snr_list, cmap='seismic',
                edgecolor='grey', linewidth=0.05)
        im.set_clim(min(snr_list), max(snr_list))
        ax.set_xlim([min(cc_dt_list), max(cc_dt_list)])
        ax.set_ylim([min(cc_max_list), max(cc_max_list)])
        ax.set_xlabel("cc_dt")
        ax.set_ylabel("cc_max")
        #add colorbar
        cbar_ax = fig.add_axes([0.95, 0.1, 0.005, 0.2])
        fig.colorbar(im, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=9)
        cbar_ax.set_xlabel('SNR(dB)', fontsize=9)
        cbar_ax.xaxis.set_label_position('top')

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
            out_file = '%s.%s.pdf' % (event_id, window_id)
        fig.savefig(out_file, format='pdf')
        #fig.savefig("misfit.pdf", bbox_inches='tight', format='pdf')


    def plot_seismograms_for_one_event(self, event_id, 
            azbin=10, win=[0,500], rayp=10,
            obs_dir='obs', syn_dir='syn', syn_band_code='MX',
            syn_suffix='.sem.sac', filter_padlen=100,
            use_prefilt=True, prefilt=('cosine', [0.11, 0.09]),
            use_STF=False, STF=('triangle', 0.0) ):
        """ Plot seismograms for one event
            azbin:
                azimuthal bin size
        """
        event = self.data['events'][event_id]
        cmt = event['gcmt']
        centroid_time = UTCDateTime(cmt['centroid_time'])

        #====== calculate traveltime curve
        model = TauPyModel(model="ak135")
        dist_ttcurve = np.arange(0.0,40,0.1)
        ttcurve_p = []
        ttcurve_P = []
        ttcurve_s = []
        ttcurve_S = []
        for dist in dist_ttcurve:
            arrivals = model.get_travel_times(source_depth_in_km=cmt['depth'],
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
        
        #====== plot seismograms in azimuthal bins
        stations = self.data['events'][event_id]['stations']
        for az in np.arange(0, 360, azbin):
            azmin = az
            azmax = az + azbin
        
            # read available stations
            data = {}
            for station_id in stations:
                station = stations[station_id]
                meta = station['meta']
                azimuth = meta['azimuth']
                if azimuth < azmin or azimuth > azmax:
                    continue
                try:
                    startime, times, detla, syn_ENZ, obs_ENZ, _, _ = \
                        self.read_seismograms_for_one_station(event_id,
                                station_id, obs_dir=obs_dir, syn_dir=syn_dir,
                                syn_band_code=syn_band_code,
                                syn_suffix=syn_suffix,
                                filter_padlen=filter_padlen, 
                                use_prefilt=use_prefilt, prefilt=prefilt, 
                                use_STF=use_STF, STF=STF)
                except:
                    continue
                # project obs/syn_ENZ to RTZ
                syn_RTZ = syn_ENZ
                obs_RTZ = obs_ENZ
                Raz = (meta['back_azimuth'] + 180.0) % 360.0
                sin_Raz = np.sin(np.deg2rad(Raz))
                cos_Raz = np.cos(np.deg2rad(Raz))
                proj_matrix = [ [ sin_Raz, cos_Raz],
                                [-cos_Raz, sin_Raz] ]
                syn_RTZ[0:2,:] = np.dot(proj_matrix, syn_ENZ[0:2,:])
                obs_RTZ[0:2,:] = np.dot(proj_matrix, obs_ENZ[0:2,:])
                # record results
                data[station_id] = {'meta':meta, 'starttime':startime, 
                        'times':times, 'delta':detla, 
                        'syn':syn_RTZ, 'obs':obs_RTZ}

            if not data:
                continue

            # create figure and axes
            fig = plt.figure(figsize=(11,8.5))
            str_title = 'azimuth(deg): {:.1f}, {:.1f}'.format(azmin, azmax)
            fig.text(0.5, 0.95, str_title, size='x-large', horizontalalignment='center')
            ax_RTZ = []
            for i in range(3):
                ax_origin = [0.05+0.3*i, 0.1]
                ax_size = [0.25, 0.8]
                ax_RTZ.append(fig.add_axes(ax_origin + ax_size))

            # make figure
            y = [ x['meta']['dist_degree'] for x in data.itervalues() ]
            ny = len(y)
            ymax = max(y) + 0.5
            ymin = min(y) - 0.5
            dy = 0.5*(ymax-ymin)/ny
        
            for sta_id in data:
                sta = data[sta_id]
                t = sta['times']
                dist_degree = sta['meta']['dist_degree']
                reduced_time = dist_degree * rayp
        
                t0 =  centroid_time - sta['starttime']
                plot_t0 = t0 + win[0] + reduced_time
                plot_t1 = t0 + win[1] + reduced_time
                idx = (t > plot_t0) & (t < plot_t1)
                    
                t_plot = t[idx] - reduced_time
                obs_RTZ = sta['obs']
                syn_RTZ = sta['syn']
                cmp_names = ['R', 'T', 'Z']
                for i in range(3):
                    #normalize data
                    obs = obs_RTZ[i, idx]
                    obs = dy*obs/max(abs(obs))
                    syn = syn_RTZ[i, idx]
                    syn = dy*syn/max(abs(syn))

                    ax = ax_RTZ[i]
                    ax.plot(t_plot, obs+dist_degree, 'k-', linewidth=0.5)
                    ax.plot(t_plot, syn+dist_degree, 'r-', linewidth=0.5)
        
                    #plot ttcurve
                    ax.plot([x[1]-rayp*x[0] for x in ttcurve_p], [x[0] for x in ttcurve_p], 'b-', linewidth=0.2)
                    ax.plot([x[1]-rayp*x[0] for x in ttcurve_P], [x[0] for x in ttcurve_P], 'b-', linewidth=0.2)
                    ax.plot([x[1]-rayp*x[0] for x in ttcurve_s], [x[0] for x in ttcurve_s], 'c-', linewidth=0.2)
                    ax.plot([x[1]-rayp*x[0] for x in ttcurve_S], [x[0] for x in ttcurve_S], 'c-', linewidth=0.2)
        
                    # config
                    ax.set_xlim(win[0], win[1])
                    ax.set_ylim(ymin, ymax)
                    ax.set_title(cmp_names[i])
                    ax.set_xlabel('t - {:.1f}*dist (s)'.format(rayp))
        
                    # ylabel 
                    if i == 0:
                        ax.set_ylabel('dist (deg)')
                    else:
                        ax.set_yticklabels([])
                    #annotate station names 
                    if i == 2:
                        ax.text(win[1], dist_degree, ' '+sta_id, \
                                verticalalignment='center', fontsize=9)
                #for i in range(3):
            #for sta_id in data:
        
            plt.show()
            #plt.savefig(fig, )
        
        # END