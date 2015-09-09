#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Compare waveform in the specific time windows

Usage:
    measure_window.py [-w WINDOW_FILE] [--obs=OBS_DIR] [--syn=SYN_DIR]
    [--meta=STATION_META] [--hdur=HALF_DURATION] [-o OUTPUT_FILE]

"""
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack, signal
from math import sin, cos
from obspy import UTCDateTime, read
import obspy.signal

#====== parameters

#window_file = 'window_test.list'
#obs_dir = 'C201002180113A_obs/sac_dis'
#syn_dir = 'C201002180113A_syn'
#
## station metadata
##station_metadata_file = 'IRIS.station'
#station_metadata_file = 'IRIS.station.correctYP'
#
## SEM channel code
#sem_channel = 'BX'
#
## filter
#freqmin = 0.015
#freqmax = 0.1
#ncorners = 2
#
## source wavelet: symmetry triangle
#half_duration = 1 # seconds
#
## window taper
#taper_width = 0.5
#
## output file
#result_file = 'result.list'

window_file           = str(sys.argv[1])
obs_dir               = str(sys.argv[2])
syn_dir               = str(sys.argv[3])
station_metadata_file = str(sys.argv[4])
sem_channel_code      = str(sys.argv[5])
freqmin               = float(sys.argv[6])
freqmax               = float(sys.argv[7])
ncorners              = int(sys.argv[8])
half_duration         = float(sys.argv[9])
taper_width           = float(sys.argv[10])
result_file           = str(sys.argv[11])

#====== read window list
with open(window_file, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.replace('\n','').split('|') for x in lines]
window_dicts = [ 
        {'network':     x[0],
         'station':     x[1],
         'location':    x[2],
         'channel':     x[3],
         'azimuth':     float(x[4]),
         'dip':         float(x[5]),
         'starttime':   UTCDateTime(x[6]),
         'endtime':     UTCDateTime(x[7]) } for x in lines ]

#TODO gather channel/az/dip/start/end to the same net/sta/loc
#   this helps in creating adjoint source at each station.

#====== read station metadata
with open(station_metadata_file, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.replace('\n','').split('|') for x in lines]
station_dicts = [
    { 'network':     x[0], 
      'station':     x[1],
      'location':    x[2],
      'channel':     x[3],
      'latitude':    float(x[4]),
      'longitude':   float(x[5]),
      'elevation':   float(x[6]),
      'depth':       float(x[7]),
      'azimuth':     float(x[8]),
      'dip':         float(x[9]),
      'starttime':   UTCDateTime(x[15]),
      'endtime':     UTCDateTime(x[16]) } for x in lines ]

#====== measure each time window
result_dicts = []
for win in window_dicts:

    # file basename (without the orietation code) 
    obs_fnbase = '{:s}.{:s}.{:s}.{:2s}'.format(win['network'], 
            win['station'], win['location'], win['channel'][0:2])
    syn_fnbase = '{:s}.{:s}.{:s}.{:2s}'.format(win['network'], 
            win['station'], win['location'], sem_channel_code[0:2])

    # check if all 3 component files exist for observed seismograms
    obs_fn = []
    obs_chan = []
    if os.path.isfile(obs_dir+'/'+obs_fnbase+'Z'):
        obs_fn.append(obs_fnbase+'Z')
        obs_chan.append(win['channel'][0:2]+'Z')
    if (os.path.isfile(obs_dir+'/'+obs_fnbase+'N') and 
            os.path.isfile(obs_dir+'/'+obs_fnbase+'E')):
        obs_fn.append(obs_fnbase+'N')
        obs_fn.append(obs_fnbase+'E')
        obs_chan.append(win['channel'][0:2]+'N')
        obs_chan.append(win['channel'][0:2]+'E')
    elif (os.path.isfile(obs_dir+'/'+obs_fnbase+'1') and 
            os.path.isfile(obs_fnbase+'2')):
        obs_fn.append(obs_fnbase+'1')
        obs_fn.append(obs_fnbase+'2')
        obs_chan.append(win['channel'][0:2]+'1')
        obs_chan.append(win['channel'][0:2]+'2')

    if len(obs_fn) != 3:
        print '[WARN] 3 components not found for obs', obs_fnbase 
        continue

    # get metadata for each observation channel
    obs_meta = []
    for i in range(3):
        channel =  [x for x in station_dicts if ( 
            x['network']  == win['network'] and 
            x['station']  == win['station'] and 
            x['location'] == win['location'] and 
            x['channel']  == obs_chan[i] and 
            x['starttime'] <= win['starttime'] and 
            x['endtime'] >= win['endtime']) ]
        if len(channel) == 1:
            obs_meta.append(channel[0])
            # check if all the channels have the same location: lat/lon/ele/dep
            if i == 0:
                lat0 = obs_meta[0]['latitude']
                lon0 = obs_meta[0]['longitude']
                ele0 = obs_meta[0]['elevation']
                dep0 = obs_meta[0]['depth']
            if i > 0:
                lat1 = obs_meta[i]['latitude']
                lon1 = obs_meta[i]['longitude']
                ele1 = obs_meta[i]['elevation']
                dep1 = obs_meta[i]['depth']
                if (lat1 != lat0 or lon1 != lon0 or ele1 != ele0 
                        or dep1 != dep0):
                    print '[WARN] different channel location for ', obs_fn, win
                    print '[WARN]   old location ', lat0,lon0,ele0,dep0
                    print '[WARN]   new location ', lat1,lon1,ele1,dep1
                    print '[WARN] Skip this time window ', win
                    continue 
        else:
            print '[WARN] failed to determine one unique metadata for ', \
                    obs_fn[i], win 
            print '[WARN] Skip this time window ', win
            continue

    # check if all 3 component files exist for synthetic seismograms
    syn_fn = []
    if (os.path.isfile(syn_dir+'/'+syn_fnbase+'Z') 
            and os.path.isfile(syn_dir+'/'+syn_fnbase+'N') 
            and os.path.isfile(syn_dir+'/'+syn_fnbase+'E')):
        syn_fn.append(syn_fnbase+'Z')
        syn_fn.append(syn_fnbase+'N')
        syn_fn.append(syn_fnbase+'E')

    if len(syn_fn) != 3:
        print '[WARN] 3 components not found for syn', syn_fnbase 
        continue

    # read observed/synthetic seismograms
    obs_seis  = read(obs_dir+'/'+obs_fn[0])
    obs_seis += read(obs_dir+'/'+obs_fn[1])
    obs_seis += read(obs_dir+'/'+obs_fn[2])

    syn_seis  = read(syn_dir+'/'+syn_fn[0])
    syn_seis += read(syn_dir+'/'+syn_fn[1])
    syn_seis += read(syn_dir+'/'+syn_fn[2])

    # get basic properties from the synthetics
    sampling_rate = syn_seis[0].stats.sampling_rate
    dt = syn_seis[0].stats.delta
    starttime = syn_seis[0].stats.starttime
    endtime = syn_seis[0].stats.endtime
    npts = syn_seis[0].stats.npts

    # make the source time function
    # NOTE now width only precise to sampling interval
    #   change to frequency domain may be better
    nstf = int(2*half_duration*sampling_rate)
    if nstf <= 1: # when window width smaller than sampling interval
        nstf = 1
    stf = signal.get_window('triang', int(2*half_duration*sampling_rate))
    stf /= sum(stf)

    # detrend(obs)/resample(obs)/shift(obs)/filter seismograms
    # note: no need to detrend synthetics, sometimes this would introduce
    # large artefacts at the begin/end of the trace
    for i in range(3):
        obs_seis[i].trim(starttime,endtime,pad=False,nearest_sample=False)
        obs_seis[i].detrend('linear')
        obs_seis[i].taper(max_percentage=0.05)
        obs_seis[i].trim(starttime,endtime,pad=True,fill_value=0.0,
                nearest_sample=False)
        obs_seis[i].resample(sampling_rate) # spectral interpolation
        # shift obs to match time samples in synthetics exactly 
        tshift = starttime - obs_seis[i].stats.starttime
        tlen = obs_seis[i].stats.npts * obs_seis[i].stats.delta
        obs_seis[i].data = fftpack.shift(obs_seis[i].data, tshift, tlen)
        obs_seis[i].stats.starttime = starttime
        obs_seis[i].trim(starttime,endtime,pad=True,fill_value=0.0,
                nearest_sample=True)
        # convolve source term on synthetics
        # NOTE: mode='same' aligns the peak/center of triangle at the beginning
        #   of the data points for the first point in the output
        syn_seis[i].data = signal.convolve(syn_seis[i].data, stf, mode='same')
        # filter seismograms
        obs_seis[i].filter('bandpass',freqmin=freqmin,freqmax=freqmax,
                corners=ncorners,zerophase=True)
        syn_seis[i].filter('bandpass',freqmin=freqmin,freqmax=freqmax,
                corners=ncorners,zerophase=True)

    # rotate to ZNE (obs)
    # TODO check the rotation formula
    obs_ZNE = np.zeros((3, npts))
    obs_ZNE[0,:], obs_ZNE[1,:], obs_ZNE[2,:] = obspy.signal.rotate.rotate2ZNE( 
            obs_seis[0].data, obs_meta[0]['azimuth'], obs_meta[0]['dip'],
            obs_seis[1].data, obs_meta[1]['azimuth'], obs_meta[1]['dip'],
            obs_seis[2].data, obs_meta[2]['azimuth'], obs_meta[2]['dip'])

    syn_ZNE = np.zeros((3, npts))
    for i in range(3):
        syn_ZNE[i,:] = syn_seis[i].data

    # make window selection function (projection, taper)
    b = win['starttime'] - starttime
    e = win['endtime'] - starttime
    ib = int(b * sampling_rate)
    ie = int(e * sampling_rate) + 1
    if ib < 0:
        ib = 0
    if ie > npts:
        ie = npts
    
    win_ZNE = np.zeros((3, npts))
    win_ZNE[:,ib:ie] = signal.tukey(ie-ib, alpha=taper_width, sym=True)

    if win['channel'][2] in ['Z', 'R', 'T']:
        sin_az = sin(np.deg2rad(win['azimuth']))
        cos_az = cos(np.deg2rad(win['azimuth']))
        sin_dip = sin(np.deg2rad(win['dip']))
        cos_dip = cos(np.deg2rad(win['dip']))
        win_ZNE[0,:] *= - sin_dip
        win_ZNE[1,:] *= cos_dip * cos_az 
        win_ZNE[2,:] *= cos_dip * sin_az 
    elif win['channel'][2] == 'H':
        win_ZNE[0,:] = 0.0
    else:
        if win['channel'][2] != 'A':
            print '[ERROR] unrecognized orientation code: ', win['channel']
            sys.exit()

    # apply window function
    obs_ZNE *= win_ZNE
    syn_ZNE *= win_ZNE

    # measure misfit between observed and synthetic seismograms
    # cross-correlation
    norm2_obs = np.sum(obs_ZNE*obs_ZNE)
    norm_obs = np.sqrt(norm2_obs)
    norm_syn = np.sqrt(np.sum(syn_ZNE*syn_ZNE))
    cc = np.zeros(npts)
    for i in range(3):
        cc += np.correlate(obs_ZNE[i,:], syn_ZNE[i,:], 'same')
    cc /= norm_obs*norm_syn
    # tshift>0: synthetic is shifted to the positive time direction
    tshift = (- npts//2 + np.arange(npts)) * dt
    # time shift at the maximum correlation
    imax = np.argmax(cc)
    cc_max = cc[imax]
    tshift_max = tshift[imax]
    # zero-lag
    cc_0 = cc[npts//2]
    ar_0 = cc_0*norm_syn/norm_obs # amplitude ratio syn/obs 

    # store result
    result = win
    result['tshift'] = tshift_max
    result['cc_max'] = cc_max
    result['cc_0'] = cc_0
    result['ar_0'] = ar_0
    result_dicts.append(result)

    print result

    # plot for debug
#   t = syn_seis[0].times()
#   plt.figure(1)
#   for i in range(3):
#       plt.subplot(411+i)
#       if i == 0:
#           plt.title(obs_fn[0])
#       plt.plot(t, obs_ZNE[i,:], 'k-', t, syn_ZNE[i,:], 'r-')
#       plt.xlim((min(t), max(t)))

#   plt.subplot(414)
#   plt.plot(tshift, cc, 'k-')
#   plt.xlim((min(tshift), max(tshift)))
#   plt.legend(['dt {:f} cc_max {:f} cc0 {:f}'.format(tshift_max, cc_max,
#       cc_0)])

#   plt.show()

#======  output results
with open(result_file, 'w') as f:
    f.write('# network|station|location|channel|azimuth|dip|'\
            'starttime|endtime|tshift_ccmax|cc_max|cc_zerolag|amplitude_ratio\n')
    for res in result_dicts:
        f.write('{:s}|{:s}|{:s}|{:s}|{:.1f}|{:.1f}|{:s}|{:s}|'\
                '{:.3f}|{:.3f}|{:.3f}|{:.3f}\n'.format(
                    res['network'], res['station'], res['location'],
                    res['channel'], res['azimuth'], res['dip'],
                    res['starttime'], res['endtime'],
                    res['tshift'], res['cc_max'], res['cc_0'], res['ar_0']))

#   print '{:s}|{:s}|{:s}|{:s}|{:.1f}|{:.1f}|{:s}|{:s}|'\
#         '{:.3f}|{:.3f}|{:.3f}|{:.3f}'.format(
#           res['network'], res['station'], res['location'], res['channel'],
#           res['azimuth'], res['dip'], res['starttime'], res['endtime'],
#           res['tshift'], res['cc_max'], res['cc_0'], res['ar_0'])
