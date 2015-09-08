#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
compare waveform in the specific time windows

"""
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack, signal
from math import sin, cos
from obspy import UTCDateTime, read
import obspy.signal

window_file = 'window.list'
obs_dir = 'C201002180113A_obs/sac_dis'
syn_dir = 'C201002180113A_syn'

# station metadata
station_metadata_file = 'IRIS.station'

# SEM
sem_channel = 'BX'
# filter
freqmin = 0.015
freqmax = 0.1
ncorners = 2

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
for win in window_dicts:

    # file basename (without the orietation code) 
    obs_fnbase = '{:s}.{:s}.{:s}.{:2s}'.format(win['network'], 
            win['station'], win['location'], win['channel'][0:2])
    syn_fnbase = '{:s}.{:s}.{:s}.{:2s}'.format(win['network'], 
            win['station'], win['location'], sem_channel[0:2])

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

    # get metadata of each observation channel
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
        else:
            print '[ERROR] in determining metadata for ', \
                    obs_fn, win 
            sys.exit()

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
        # filter seismograms
        obs_seis[i].filter('bandpass',freqmin=freqmin,freqmax=freqmax,
                corners=ncorners,zerophase=True)
        syn_seis[i].filter('bandpass',freqmin=freqmin,freqmax=freqmax,
                corners=ncorners,zerophase=True)

    # rotate to ZNE (obs)
    print obs_meta
    sys.exit()
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
    win_ZNE[:,ib:ie] = signal.tukey(ie-ib, alpha=0.5, sym=True)

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
    #obs_ZNE *= win_ZNE
    #syn_ZNE *= win_ZNE

    t = syn_seis[0].times()

    plt.figure(1)
    for i in range(3):
        plt.subplot(311+i)
        plt.plot(t, obs_ZNE[i,:], 'k-', t, syn_ZNE[i,:], 'r-')
    plt.show()
