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
import scipy.signal #fftconvolve, tukey
import matplotlib.pyplot as plt
from obspy import UTCDateTime, read

import lanczos_interp1 as interp

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

# lanczos interpolation
na = 10 # sinc window width

# cross-correlation
cc_time_step = 0.01
cc_max_time_shift = 10

# read from command line
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

fp_result = open(result_file,'w')

#====== read window list

#TODO put comment lines into the result file
#   this will keep the event info and window specfications
with open(window_file, 'r') as f:
    for x in f.readlines():
        if x.startswith('#'):
            fp_result.write(x)

with open(window_file, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.replace('\n','').split('|') for x in lines]
windows = {}
for x in lines:
    net_sta_loc_comp3 = (x[0], x[1], x[2], x[3], x[4], x[5],
            x[6], x[7], x[8], x[9])
    win = {'phase':       x[10],
           'orientation': x[11],
           'azimuth':     float(x[12]),
           'dip':         float(x[13]),
           'starttime':   UTCDateTime(x[14]),
           'endtime':     UTCDateTime(x[15])}
    if net_sta_loc_comp3 not in windows:
        windows[net_sta_loc_comp3] = []
    windows[net_sta_loc_comp3].append(win)

#====== measure each time window
for net_sta_loc_comp3 in windows:

    net_sta_loc = '.'.join(net_sta_loc_comp3[0:3])
    obs_comp3 = net_sta_loc_comp3[3].split(',')
    obs_az3 = [float(x) for x in net_sta_loc_comp3[4].split(',')]
    obs_dip3 = [float(x) for x in net_sta_loc_comp3[5].split(',')]

    # file names: obs, syn
    obs_fn = ['{:s}/{:s}.{:s}'.format(obs_dir, net_sta_loc, x) 
            for x in obs_comp3]
    syn_fn = ['{:s}/{:s}.{:s}{:s}'.format(
        syn_dir, net_sta_loc, sem_channel_code[0:2], x) 
        for x in ['Z', 'N', 'E']]

    # check file existence: obs, syn
    flag_skip = False
    for i in range(3):
        if not os.path.isfile(obs_fn[i]):
            print '[WARN] file not found: ', obs_fn[i]
            flag_skip = True
            break
        if not os.path.isfile(syn_fn[i]):
            print '[WARN] file not found: ', syn_fn[i]
            flag_skip = True
            break 
    if flag_skip:
        print '[WARN] Skip ', net_sta_loc_comp3
        continue

    # read observed/synthetic seismograms
    obs_seis  = read(obs_fn[0])
    obs_seis += read(obs_fn[1])
    obs_seis += read(obs_fn[2])

    syn_seis  = read(syn_fn[0])
    syn_seis += read(syn_fn[1])
    syn_seis += read(syn_fn[2])

    # get basic properties from the synthetics
    # assume the same for all 3-components of syn
    # TODO check 3-comp of syn have same npts/dt/starttime/endtime
    sampling_rate = syn_seis[0].stats.sampling_rate
    dt = syn_seis[0].stats.delta
    starttime = syn_seis[0].stats.starttime
    endtime = syn_seis[0].stats.endtime
    npts = syn_seis[0].stats.npts
    t = np.arange(npts)*dt

    # source time function
    freq = np.fft.rfftfreq(npts, d=dt)
    # triangle(two equal sides)
    src_spectrum = np.sinc(freq*half_duration)**2

    # detrend(obs)/resample(obs)/shift(obs)/filter seismograms
    # note: no need to detrend synthetics, sometimes this would introduce
    # large artefacts at the begin/end of the trace
    for i in range(3):
        # filter obs
        obs_seis[i].detrend('linear')
        obs_seis[i].taper(max_percentage=0.05)
        obs_seis[i].filter('bandpass',freqmin=freqmin,freqmax=freqmax,
                corners=ncorners,zerophase=True)
        # interpolate obs to the same time samples as syn
        obs_seis[i].data = interp.lanczos_interp1(
                obs_seis[i].data, obs_seis[i].stats.delta,
                t+(starttime-obs_seis[i].stats.starttime), na)
        # convolve synthetics with source time function 
        syn_seis[i].data = np.fft.irfft(
                np.fft.rfft(syn_seis[i].data)*src_spectrum)
        # filter seismograms
        syn_seis[i].filter('bandpass',freqmin=freqmin,freqmax=freqmax,
                corners=ncorners,zerophase=True)

    # rotate 3-comp of obs to ZNE coordinates
    obs_ZNE = np.zeros((3, npts))
    # projection matrix: obs = proj * ZNE
    proj_matrix = np.zeros((3, 3))
    for i in range(3):
        sin_az = np.sin(np.deg2rad(obs_az3[i]))
        cos_az = np.cos(np.deg2rad(obs_az3[i]))
        sin_dip = np.sin(np.deg2rad(obs_dip3[i]))
        cos_dip = np.cos(np.deg2rad(obs_dip3[i]))
        # column vector = obs orientation
        proj_matrix[0,i] = -sin_dip
        proj_matrix[1,i] = cos_dip * cos_az
        proj_matrix[2,i] = cos_dip * sin_az
        obs_ZNE[i,:] = obs_seis[i].data
    # inverse projection matrix: ZNE = inv(proj) * obs
    obs_ZNE = np.dot(np.linalg.inv(proj_matrix), obs_ZNE)

    # synthetic data matrix
    syn_ZNE = np.zeros((3, npts))
    for i in range(3):
        syn_ZNE[i,:] = syn_seis[i].data

    # measure each time window
    obs_ZNE_win = np.zeros((3, npts))
    syn_ZNE_win = np.zeros((3, npts))
    taper = np.zeros((3, npts))

    for win in windows[net_sta_loc_comp3]:
        # make taper window
        tb = win['starttime'] - starttime
        te = win['endtime'] - starttime
        ib = int(tb * sampling_rate)
        ie = int(te * sampling_rate) + 1
        if ib < 0: ib = 0
        if ie > npts: ie = npts
        taper[:,:] = 0.0 #NOTE don't forget to reset all the values to zero
        taper[:,ib:ie] = scipy.signal.tukey(ie-ib, alpha=taper_width, sym=True)
    
        # projection matrix
        proj_matrix[:,:] = 0.0 #reselt to zero
        if win['orientation'] in ['Z', 'R', 'T']:
            sin_az = np.sin(np.deg2rad(win['azimuth']))
            cos_az = np.cos(np.deg2rad(win['azimuth']))
            sin_dip = np.sin(np.deg2rad(win['dip']))
            cos_dip = np.cos(np.deg2rad(win['dip']))
            n = np.array([[-sin_dip],[cos_dip * cos_az],[cos_dip * sin_az]])
            proj_matrix = np.dot(n, n.transpose())
        elif win['orientation'] == 'H': # horizontal vector 2d
            proj_matrix[0,0] = 0.0
            proj_matrix[1,1] = 1.0
            proj_matrix[2,2] = 1.0
        elif win['orientation'] == 'F': # full 3d vector
            proj_matrix[0,0] = 1.0
            proj_matrix[1,1] = 1.0
            proj_matrix[2,2] = 1.0
        else:
            print '[WARN] unrecognized orientation code: ', \
                    win['orientation'], net_sta_loc
            break

        # apply taper and projection
        obs_ZNE_win = np.dot(proj_matrix, obs_ZNE) * taper
        syn_ZNE_win = np.dot(proj_matrix, syn_ZNE) * taper
    
        # measure misfit between observed and synthetic seismograms
        # cross-correlation
        norm2_obs = np.sum(obs_ZNE_win**2)
        norm_obs = np.sqrt(norm2_obs)
        norm_syn = np.sqrt(np.sum(syn_ZNE_win**2))
        if norm_obs == 0:
            print '[WARN] Skip empty obs trace: ', win, obs_fn
            continue
        cc = np.zeros(2*npts-1)
        for i in range(3):
            # NOTE the order (obs,syn) is important. The positive time on 
            #   CC means shifting syn in the positive time direction
            cc += scipy.signal.fftconvolve(
                    obs_ZNE_win[i,:], syn_ZNE_win[i,::-1], 'full')
        cc /= norm_obs*norm_syn
        # tshift>0: synthetic is shifted along the positive time direction
        ncc = int(cc_max_time_shift/cc_time_step)
        tshift = np.arange(-ncc,ncc+1)*cc_time_step
        # interpolate cc to finer time samples
        if dt < cc_time_step:
            print '[WARN] syn_dt(%f) < cc_time_step(%f)\n' % (dt, cc_time_step)
        ti = (npts-1)*dt + tshift # -(npts-1)*dt: begin time in cc
        cci = interp.lanczos_interp1(cc, dt, ti, na)
        # time shift at the maximum correlation
        imax = np.argmax(cci)
        cc_max = cci[imax]
        tshift_cc_max = tshift[imax]
        # zero-lag correlation
        cc_0 = cci[ncc]
        ar_0 = cc_0*norm_syn/norm_obs # amplitude ratio syn/obs 
    
        # output results
        fp_result.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|' \
                '%s|%s|%.1f|%.1f|%s|%s|%.3f|%.3f|%.3f|%.3f\n' % (
                    net_sta_loc_comp3[0], net_sta_loc_comp3[1], 
                    net_sta_loc_comp3[2], net_sta_loc_comp3[3],
                    net_sta_loc_comp3[4], net_sta_loc_comp3[5],
                    net_sta_loc_comp3[6], net_sta_loc_comp3[7],
                    net_sta_loc_comp3[8], net_sta_loc_comp3[9],
                    win['phase'], win['orientation'], 
                    win['azimuth'], win['dip'],
                    win['starttime'], win['endtime'],
                    tshift_cc_max, cc_max, cc_0, ar_0))
 
        # DEBUG 
#       t = syn_seis[0].times()
#       for i in range(3):
#           plt.subplot(411+i)
#           if i == 0:
#               plt.title(net_sta_loc)
#           plt.plot(t, obs_ZNE_win[i,:], 'k-', t, syn_ZNE_win[i,:], 'r-')
#           plt.xlim((min(t), max(t)))
#           plt.ylabel(obs_comp3[i])
#   
#       plt.subplot(414)
#       plt.plot(tshift, cci, 'k-')
#       plt.xlim((min(tshift), max(tshift)))
#       plt.ylabel('{:s} {:s}'.format(win['phase'], win['orientation']))
#       plt.legend(['dt {:.2f} cc_max {:.3f} cc0 {:.3f} ar0 {:.3f}'.format(
#           tshift_cc_max, cc_max, cc_0, ar_0)])
#   
#       plt.show()

    # end for win in windows[net_sta_loc_comp3]
#   sys.exit()

# end for net_sta_loc_comp3 in windows

fp_result.close()

# END 
