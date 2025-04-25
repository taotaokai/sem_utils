#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import json
import numpy as np

# read command line args
data_dir = str(sys.argv[1])
obs_dir = str(sys.argv[2])
syn_dir = str(sys.argv[3])
freqmin = float(sys.argv[4])
freqmax = float(sys.argv[5])

#------ input files
CMTSOLUTION_file = '%s/CMTSOLUTION.obs' % (data_dir)
station_file = '%s/station.txt' % (data_dir)

#------
print "\ninitialize\n"
misfit = Misfit()

#------
print "\n====== setup event\n"
misfit.setup_event_from_CMTSOLUTION(CMTSOLUTION_file)
event_id = [ key for key in misfit.data['events'] ][0]

#------
print "\n====== setup stations\n"
misfit.setup_stations_from_metafile(station_file)

#------
print "\n====== setup windows\n"
window_list = [ 
        ('F','p,P',[-30,50]), 
        ('F','s,S',[-40,70]) ]
print "window_list= ", window_list
filter_param=('butter', 3, [freqmin, freqmax])
print "filter_param= ", filter_param
misfit.setup_windows(window_list=window_list,
        filter_param=filter_param)

#------
print "\n====== measure windows\n"
#window_id_list = ['F.p,P','F.s,S']
window_id_list = [ '.'.join(x[0:2]) for x in window_list ]
weight_param={'SNR':[10, 15], 'CCmax':[0.6,0.8], 'CC0':[0.5,0.7]}
print "weight_param= ", weight_param

event = misfit.data['events'][event_id]
gcmt = event['gcmt']

hdur_list = np.arange(0, 16)
wCCmax_list = np.zeros(16)

for i in range(16):
    gcmt['half_duration'] = hdur_list[i]
    misfit.measure_windows_for_one_event(event_id=event_id,
        obs_dir=obs_dir, syn_dir=syn_dir,
        use_STF=True, plot=False, output_adj=False,
        weight_param=weight_param)
    wCCmax = 0.0
    weight = 0.0
    stations = event['stations']
    for station_id in stations:
        station = stations[station_id]
        if station['stat']['code'] < 0:
            continue
        windows = station['windows']
        for window_id in windows:
            window = windows[window_id]
            if window['stat']['code'] != 1:
                continue
            weight += window['weight']
            wCCmax += window['weight'] * window['misfit']['CCmax']

    wCCmax_list[i] = wCCmax / weight
    

print "#====== half-duration search results\n"
for i in range(16):
    print hdur_list[i], wCCmax_list[i]