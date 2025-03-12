#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import warnings
import re
#
import numpy as np
import scipy.signal as signal
import tables
#
from obspy import UTCDateTime, read, Trace, geodetics
from obspy.taup import TauPyModel
from obspy.imaging.beachball import beach
#
from matplotlib import colors, ticker, cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def spheredist(lat0,lon0,lat1,lon1):
    d2r = np.pi/180.0
    lat0, lon0 = lat0*d2r, lon0*d2r
    lat1, lon1 = lat1*d2r, lon1*d2r
    # v0 = [np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
    # v1 = [np.cos(lat1)*np.cos(lon1), np.cos(lat1)*np.sin(lon1), np.sin(lat1)]
    return np.arccos(np.sin(lat0)*np.sin(lat1) + np.cos(lat0)*np.cos(lat1)*np.cos(lon1-lon0))

def centerMap(lats,lons,scale):
    """
    https://stackoverflow.com/questions/13240642/automatically-center-matplotlib-basemap-onto-data
    """
    # NOTE: only works for northern hemisphere!
    #Assumes -90 < Lat < 90 and -180 < Lon < 180, and
    # latitude and logitude are in decimal degrees
    earthRadius = 6378100.0 #earth's radius in meters
    northLat = max(lats)
    southLat = min(lats)
    westLon = max(lons)
    eastLon = min(lons)
    # average between max and min longitude
    # lon0 = ((westLon-eastLon)/2.0)+eastLon
    lon0 = np.median(lons)
    # a = the height of the map
    b = spheredist(northLat,westLon,northLat,eastLon)*earthRadius/2
    c = spheredist(northLat,westLon,southLat,lon0)*earthRadius
    # use pythagorean theorom to determine height of plot
    mapH = pow(pow(c,2)-pow(b,2),1./2)
    arcCenter = (mapH/2)/earthRadius
    lat0 = southLat + arcCenter*180./np.pi
    # distance between max E and W longitude at most souther latitude
    # widest part on the map, either at the equator or at the latitude closest to the equator
    minLat = min(abs(southLat), abs(northLat))
    if np.sign(southLat) != np.sign(northLat): minLat = 0
    mapW = spheredist(minLat,westLon,minLat,eastLon)*earthRadius
    return lat0,lon0,mapW*scale,mapH*scale


#
def plot_seismogram_3comp(event_h5grp,
                          savefig=False,
                          out_dir='plot',
                          data_tag='DATA_DISP',
                          syn_tag='SYN_DISP',
                          # azimuthal_bin=10, # bin size of azimuth (degree)
                          max_bin_az=10, # max azimuthal size of each bin
                          max_bin_sta=20, # max num of stations in each bin
                          reduce_rayp=10., # reduced time = time_after_origin - dist_degree * rayp (sec/degree)
                          time_limit=[-100, 1000], # in reduced time (sec)
                          freq_limit=[0.01, 0.1], # filter frequency band (Hz)
                          dist_limit=None, # limit epi-distance (degree)
                          min_SNR=10, # minimum signal-to-noise ratio: 20*log10(Asig/Anoise)
                          noise_twin=[-500,-100], # relative to first arrival
                          ):

    event = event_h5grp._v_attrs['event']
    event_name = event.event_descriptions[0].text

    for origin in event.origins:
        if origin.resource_id == event.preferred_origin_id:
            event_origin = origin
            break
    # print(event_origin)
    evla = event_origin.latitude
    evlo = event_origin.longitude
    evdp = event_origin.depth / 1000.0 # km

    for focmec in event.focal_mechanisms:
        if focmec.resource_id == event.preferred_focal_mechanism_id:
            event_focmec = focmec
            break
    # print(event_focmec)
    mt = event_focmec.moment_tensor.tensor
    focmec = [mt.m_rr, mt.m_tt, mt.m_pp, mt.m_rt, mt.m_rp, mt.m_tp]

    # check argument parameters
    # if azimuthal_bin <= 0.0:
    #     raise Exception(f"azimuthal_bin({azimuthal_bin}) must > 0.0")
    if max_bin_az <= 0.0:
         raise Exception(f"max_bin_az({max_bin_az}) must > 0.0")

    comp_name = ['R', 'T', 'Z']


    #----- find stations which have data in the plotting time window
    stations = []
    for net_sta_grp in event_h5grp:
        network = net_sta_grp._v_attrs['network']
        station = net_sta_grp._v_attrs['station']
        stnm = f'{network}.{station}'
        stla = net_sta_grp._v_attrs['latitude']
        stlo = net_sta_grp._v_attrs['longitude']
        dist, az, baz = geodetics.gps2dist_azimuth(evla, evlo, stla, stlo)
        dist_degree = geodetics.kilometer2degrees(dist/1000.0)
        if dist_limit:
            if dist_degree < min(dist_limit) or dist_degree > max(dist_limit):
                continue
        # arrivals = taup_model.get_travel_times(
        #     source_depth_in_km=evdp,
        #     distance_in_degree=dist_degree,
        #     phase_list=['p','P'])
        # first_arrival_time = event_origin.time + min([arr.time for arr in arrivals])
        # find the data trace containing the required plotting window
        # h5data = None
        reduced_time = dist_degree * reduce_rayp
        plot_t1 = event_origin.time + min(time_limit) + reduced_time
        plot_t2 = event_origin.time + max(time_limit) + reduced_time
        # for arr in net_sta_grp:
        #     tb = UTCDateTime(arr._v_attrs['starttime'])
        #     te = UTCDateTime(arr._v_attrs['endtime'])
        #     # if tb < plot_t1 and plot_t2 < te:
        #     if tb < plot_t1 and plot_t1 < te:
        #         h5data = arr
        #         break
        # if not h5data:
        #     print(f'[WARN] {stnm} [{tb}, {te}] does not cover plot_twin [{plot_t1}, {plot_t2}]')
        #     continue
        sta = {}
        if data_tag in net_sta_grp:
            sta['h5data'] = net_sta_grp[data_tag]
        else:
            continue
        if syn_tag in net_sta_grp:
            sta['h5syn'] = net_sta_grp[syn_tag]
        sta['latitude'] = stla
        sta['longitude'] = stlo
        sta['dist_degree'] = dist_degree
        sta['azimuth'] = az
        sta['back_azimuth'] = baz
        sta['name'] = stnm
        # sta['first_arrival'] = first_arrival_time
        stations.append(sta)

    if not stations: return

    #----- traveltime curve
    print(f'[INFO] calculate traveltime curve')
    taup_model = TauPyModel(model="ak135")
    max_dist = max([sta['dist_degree'] for sta in stations])
    dist_ttcurve = np.arange(0.0, max_dist+2, 1.0)
    ttcurve_p = []
    ttcurve_P = []
    ttcurve_s = []
    ttcurve_S = []
    first_arrivals = []
    for dist in dist_ttcurve:
        arrivals = taup_model.get_travel_times(
            source_depth_in_km=evdp,
            distance_in_degree=dist,
            phase_list=['p','P','s','S'])
        for arr in arrivals:
            if arr.name == 'p':
                ttcurve_p.append((arr.distance, arr.time, arr.ray_param))
            elif arr.name == 'P':
                ttcurve_P.append((arr.distance, arr.time, arr.ray_param))
            elif arr.name == 's':
                ttcurve_s.append((arr.distance, arr.time, arr.ray_param))
            elif arr.name == 'S':
                ttcurve_S.append((arr.distance, arr.time, arr.ray_param))
        if len(arrivals) > 0:
            first_arrivals.append(min([arr.time for arr in arrivals]))
        else:
            first_arrivals.append(0)
    # sort phases
    ttcurve_p = sorted(ttcurve_p, key=lambda x: x[2])
    ttcurve_P = sorted(ttcurve_P, key=lambda x: x[2])
    ttcurve_s = sorted(ttcurve_s, key=lambda x: x[2])
    ttcurve_S = sorted(ttcurve_S, key=lambda x: x[2])

    # interpolate first arrival times for each station
    to_remove = []
    for sta in stations:
        tt = np.interp(sta['dist_degree'], dist_ttcurve, first_arrivals)
        sta['first_arrival'] = tt
        # remove stations without enough record before first arrival to get noise amplitude
        h5data = sta['h5data']
        data_starttime = UTCDateTime(h5data._v_attrs['starttime'])
        data_t0 = data_starttime - event_origin.time
        noise_t0 = sta['first_arrival'] + min(noise_twin)
        if data_t0 > noise_t0:
            print(f'[WARN] not enough record length for noise, {data_t0} > {noise_t0}, skip {sta["name"]}')
            to_remove.append(sta)
    stations = [sta for sta in stations if sta not in to_remove]

    # map configuration
    stla_all = [sta['latitude'] for sta in stations]
    stlo_all = [sta['longitude'] for sta in stations]
    parallels = np.arange(-90.,90,20.)
    meridians = np.arange(0.,360,20.)
    lat0center,lon0center,mapWidth,mapHeight = centerMap([*stla_all, evla], [*stlo_all, evlo], 1.1)

    # plot stations waveforms (one figure for each azimuthal bin)
    stations = sorted(stations, key=lambda x: x['azimuth']) # sort by azimuth
    ista, nsta = 0, len(stations)
    num_bins = 0
    if max_bin_sta < 1: max_bin_sta = nsta
    while ista < len(stations):
        azmin = stations[ista]['azimuth']
        azmax = azmin + max_bin_az

        stations_azimuthal_bin = []
        nsta_bin = 0
        for ii in range(ista, nsta):
            sta = stations[ii]
            if sta['azimuth'] > azmax:
                ista = ii
                break
            stations_azimuthal_bin.append(sta)
            nsta_bin = nsta_bin + 1
            if nsta_bin >= max_bin_sta:
                ista = ii + 1
                break
            if ii == nsta - 1: # last station
                ista = nsta # stop while loop
        azmax = stations_azimuthal_bin[-1]['azimuth']

        print(f'[INFO] plot azimuthal bin: {azmin:05.1f} - {azmax:05.1f}')

        # skip empty azimuthal_bin
        if nsta_bin == 0:
            print(f'No station in the azimuthal bin [{azmin}, {azmax}]')
            continue

        #---- create figure
        # fig = plt.figure(figsize=(8.5, 11)) # US letter
        fig = plt.figure(figsize=(8.27, 11.69)) # A4
        str_title = f'{event_name} (az=[{azmin:05.1f}, {azmax:05.1f}], freq={freq_limit}, SNR>{min_SNR})'
        fig.text(0.5, 0.98, str_title, size='x-large', horizontalalignment='center')
        #---- plot station/event map
        ax_height = 0.25
        ax_width = ax_height * mapWidth / mapHeight
        ax_size = [ax_width, ax_height]
        ax_origin = [0.5-ax_width/2, 0.75]
        ax_map = fig.add_axes(ax_origin + ax_size)
        # ax_bm = Basemap(projection='poly', resolution='l', area_thresh=1000.,
        #     llcrnrlat=min_lat, llcrnrlon=min_lon,
        #     urcrnrlat=max_lat, urcrnrlon=max_lon,
        #     lat_0=lat_0, lon_0=lon_0, ax=ax_map)
        # ax_bm = Basemap(projection='ortho', resolution='l', lat_0=lat_0, lon_0=lon_0)
        ax_bm = Basemap(projection='stere', resolution='l', #area_thresh=10000.,
                        width=mapWidth, height=mapHeight, lat_0=lat0center, lon_0=lon0center,
                        ax=ax_map)
        ax_bm.drawcoastlines(linewidth=0.1)
        ax_bm.drawcountries(linewidth=0.1)
        ax_bm.drawlsmask()
        ax_bm.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,0], fontsize=10) #, fmt='%3.0f')
        ax_bm.drawmeridians(meridians, linewidth=0.1, labels=[0,0,0,1], fontsize=10) #, fmt='%3.0f')
        sx, sy = ax_bm(stlo_all, stla_all)
        ax_bm.scatter(sx, sy, s=5, marker='^', facecolor='blue')
        # plot focal mechanism
        sx, sy = ax_bm(evlo, evla)
        # bb_width = 110000.0 * np.abs(max(stlo_all)-min(stlo_all)) * 0.1
        bb_width = max(mapWidth, mapHeight) * 0.05
        b = beach(focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor='r')
        ax_map.add_collection(b)
        # plot the station location
        stla = [sta['latitude'] for sta in stations_azimuthal_bin]
        stlo = [sta['longitude'] for sta in stations_azimuthal_bin]
        sx, sy = ax_bm(stlo, stla)
        ax_bm.scatter(sx, sy, s=5, marker='^', facecolor='red')

        #-- create axis for seismograms
        ax_RTZ = []
        for i in range(3):
            ax_origin = [0.06+0.29*i, 0.05]
            ax_size = [0.285, 0.70]
            ax_RTZ.append(fig.add_axes(ax_origin + ax_size))

        #-- plot traveltime curves
        for i in range(3):
            ax = ax_RTZ[i]
            ax.plot([x[1]-reduce_rayp*x[0] for x in ttcurve_p], \
                [x[0] for x in ttcurve_p], 'b-', linewidth=0.1, alpha=0.2)
            ax.plot([x[1]-reduce_rayp*x[0] for x in ttcurve_P], \
                [x[0] for x in ttcurve_P], 'b-', linewidth=0.1, alpha=0.2)
            ax.plot([x[1]-reduce_rayp*x[0] for x in ttcurve_s], \
                [x[0] for x in ttcurve_s], 'c-', linewidth=0.1, alpha=0.2)
            ax.plot([x[1]-reduce_rayp*x[0] for x in ttcurve_S], \
                [x[0] for x in ttcurve_S], 'c-', linewidth=0.1, alpha=0.2)

        #-- ylim setting
        y = [sta['dist_degree'] for sta in stations_azimuthal_bin]
        ny = len(y)
        plot_dy = 0.5*(max(y)-min(y)+1)/ny
        # if dist_limit:
        #     plot_ymax = max(dist_limit) + 2*plot_dy
        #     plot_ymin = min(dist_limit) - 2*plot_dy
        # else:
        plot_ymax = max(y) + 2*plot_dy
        plot_ymin = min(y) - 2*plot_dy

        # plot each station
        for sta in stations_azimuthal_bin:
            h5data = sta['h5data']
            data_starttime = UTCDateTime(h5data._v_attrs['starttime'])
            data_sampling_rate = h5data._v_attrs['sampling_rate']
            data_npts = h5data._v_attrs['npts']
            data_orientation = h5data._v_attrs['component']
            data_endtime = data_starttime + (data_npts-1) / data_sampling_rate

            ind_zcomp = []
            ind_hcomp = []
            for i in range(len(data_orientation)):
                comp = data_orientation[i]
                if comp['code'] == b'Z':
                    ind_zcomp.append(i)
                else:
                    ind_hcomp.append(i)
            if len(ind_zcomp) != 1:
                print(f'[WARN] not exactly one vertical component: {orientation}, skip {sta}')
                continue
            if len(ind_hcomp) != 0 and len(ind_hcomp) != 2:
                print(f'[WARN] not exactly zero or two horizontal components: {orientation}, skip {sta}')
                continue

            assert(data_orientation[ind_zcomp[0]]['dip'] == -90.0)
            data_Z = h5data[ind_zcomp[0], :]
            data_RT = None
            if len(ind_hcomp) == 2:
                # rotate horizontal components to E-N
                # H12 = proj * EN => EN = inv(proj) * H12
                proj_EN_to_H12 = np.zeros((2, 2))
                i = 0
                for ind in ind_hcomp:
                    comp = data_orientation[ind]
                    az = np.deg2rad(comp['azimuth'])
                    assert(comp['dip'] == 0)
                    sin_az, cos_az = np.sin(az), np.cos(az)
                    proj_EN_to_H12[i,0] = sin_az # proj E to H
                    proj_EN_to_H12[i,1] = cos_az # proj N to H
                    i = i + 1
                # inverse projection matrix: EN = inv(proj) * H12
                proj_H12_to_EN = np.linalg.inv(proj_EN_to_H12)

                # rotate from EN to RT (T-R-Zup: right-hand convention)
                radial_az = (sta['back_azimuth'] + 180.0) % 360.0
                sin_az = np.sin(np.deg2rad(radial_az))
                cos_az = np.cos(np.deg2rad(radial_az))
                proj_EN2RT = np.array([[sin_az,  cos_az],
                                       [cos_az, -sin_az]])
                # print(proj_EN_to_H12, proj_H12_to_EN, proj_EN2RT)
                data_RT = np.dot(proj_EN2RT, np.dot(proj_H12_to_EN, h5data[ind_hcomp, :]))

            # filter
            sos = signal.butter(4, freq_limit, 'bp', fs=data_sampling_rate, output='sos')
            data_Z = signal.sosfiltfilt(sos, data_Z)
            if data_RT is not None:
                data_RT = signal.sosfiltfilt(sos, data_RT)

            # get plot time
            dist_degree = sta['dist_degree']
            reduced_time = dist_degree * reduce_rayp
            # relative time of first data sample to event origin time
            data_t0 = data_starttime - event_origin.time
            times = np.arange(data_npts) / data_sampling_rate + data_t0 # times after origin

            # noise window
            noise_t0 = sta['first_arrival'] + min(noise_twin)
            noise_t1 = sta['first_arrival'] + max(noise_twin)
            # if data_t0 > noise_t0:
            #     print(f'[WARN] not enough record length for noise, {data_t0} > {noise_t0}, skip {sta["name"]}')
            #     continue
            noise_idx = (times > noise_t0) & (times < noise_t1)

            # plot time window
            plot_t0 = min(time_limit) + reduced_time
            plot_t0 = max(plot_t0, noise_t0)
            plot_t1 = max(time_limit) + reduced_time
            plot_idx = (times > plot_t0) & (times < plot_t1)
            t_plot = times[plot_idx] - reduced_time

            # plot seismograms
            amp_Z = np.max(data_Z[plot_idx]**2)**0.5
            noise_Z = np.max(data_Z[noise_idx]**2)**0.5
            amp_RT = 0
            noise_RT = 0
            if data_RT is not None:
                amp_R = np.max(data_RT[0,plot_idx]**2)**0.5
                noise_R = np.max(data_RT[0,noise_idx]**2)**0.5
                amp_T = np.max(data_RT[1,plot_idx]**2)**0.5
                noise_T = np.max(data_RT[1,noise_idx]**2)**0.5
            # data_max_amp = (amp2_Z + amp2_RT)**0.5

            snr_Z = 20*np.log10(amp_Z/noise_Z)
            ax = ax_RTZ[2]
            h_lines = ax.plot(t_plot, plot_dy*data_Z[plot_idx]/amp_Z+dist_degree, 'k-', linewidth=0.2)
            h_text = ax.text(max(time_limit), dist_degree, f' {sta["name"]}', verticalalignment='center', fontsize=7)
            if min_SNR and snr_Z < min_SNR:
                for h in h_lines: h.set_alpha(0.1)
                h_text.set_alpha(0.2)

            if data_RT is not None:
                ax = ax_RTZ[0]
                h_lines = ax.plot(t_plot, plot_dy*data_RT[0,plot_idx]/amp_R+dist_degree, 'k-', linewidth=0.2)
                snr_R = 20*np.log10(amp_R/noise_R)
                if min_SNR and snr_R < min_SNR:
                    for h in h_lines: h.set_alpha(0.1)

                # # annotate amplitude
                # ax.text(max(time_limit), dist_degree, '%.1e ' % (amp_RT),
                #     verticalalignment='bottom',
                #     horizontalalignment='right',
                #     fontsize=7, color='black')

                ax = ax_RTZ[1]
                h_lines = ax.plot(t_plot, plot_dy*data_RT[1,plot_idx]/amp_T+dist_degree, 'k-', linewidth=0.2)
                snr_T = 20*np.log10(amp_T/noise_T)
                if min_SNR and snr_T < min_SNR:
                    for h in h_lines: h.set_alpha(0.1)

            # print(f'{sta}: snr_Z={snr_Z}, snr_R={snr_RT}')

        # set axes limits and lables, annotation
        for i in range(3):
            ax = ax_RTZ[i]
            ax.set_xlim(min(time_limit), max(time_limit))
            ax.set_ylim(plot_ymin, plot_ymax)
            ax.set_title(comp_name[i])
            ax.set_xlabel('t - {:.1f}*dist (s)'.format(reduce_rayp))
            ax.tick_params(axis='both',labelsize=10)
            # ylabel
            if i == 0:
                ax.set_ylabel('dist (deg)')
            else:
                ax.set_yticklabels([])

        # save figures
        if savefig:
            out_file = f'{out_dir}/bin{num_bins:03d}_az_{azmin:05.1f}_{azmax:05.1f}.pdf'
            plt.savefig(out_file, format='pdf')
        else:
            plt.show()
        plt.close(fig)

        num_bins = num_bins + 1

    h5f.close()


if __name__ == '__main__':
    data_file_h5 = sys.argv[1] #'syn.h5'
    figure_dir = sys.argv[2] #'figures_syn'

    if tables.is_hdf5_file(data_file_h5):
      h5f = tables.open_file(data_file_h5, mode="r")
    else:
      print(f'[Error] not hdf5 file {data_file_h5}\n')
      sys.exit()

    for evt_g in h5f.root:
        event = evt_g._v_attrs['event']
        event_name = event.event_descriptions[0].text
        out_dir = os.path.join(figure_dir, event_name)
        os.makedirs(out_dir, exist_ok=True)
        print(f'plotting {event_name}')
        plot_seismogram_3comp(evt_g, savefig=True, out_dir=out_dir,
                              data_tag='DATA_DISP',
                              max_bin_az=10, max_bin_sta=20,
                              reduce_rayp=10,
                              time_limit=[-100, 250],
                              min_SNR=10, noise_twin=[-450,-100],
                              # min_SNR=-10, noise_twin=[-10,0],
                              freq_limit=[0.02, 0.1],
                              dist_limit=[0, 40],
                              )
