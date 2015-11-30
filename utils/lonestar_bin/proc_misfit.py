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
adj_dir = str(sys.argv[4])
#adj_dir = 'adj'

#------ input files
CMTSOLUTION_file = '%s/CMTSOLUTION' % (data_dir)
station_file = '%s/station.txt' % (data_dir)

#------
print "\ninitialize\n"
misfit = Misfit()

#------
print "\nsetup event\n"
misfit.setup_event_from_CMTSOLUTION(CMTSOLUTION_file)
event_id = [ key for key in misfit.data['events'] ][0]

#------
print "\nsetup stations\n"
misfit.setup_stations_from_metafile(station_file)

#------
print "\nsetup windows\n"
window_list = [ 
        ('F','p,P',[-30,50]), 
        ('F','s,S',[-40,70]) ]
filter_param=('butter', 3, [0.012, 0.08])
misfit.setup_windows(window_list=window_list,
        filter_param=filter_param)

#------
#print "\nmeasure windows\n"
window_id_list = ['F.p,P','F.s,S']
weight_param={'SNR':[10, 15], 'cc_max':[0.6,0.8], 'cc_0':[0.5,0.7]}
#print STF
# 
misfit.measure_windows_for_one_event(event_id=event_id,
        obs_dir=obs_dir, syn_dir=syn_dir, adj_dir=adj_dir,
        adj_window_id_list=window_id_list,
        use_STF=True, plot=False, output_adj=True,
        weight_param=weight_param)

#------
print "\nsave data\n"
misfit.save(filename='%s/misfit.json' % (adj_dir))

#------
#print "\nplot seismograms\n"
#misfit.plot_seismograms(event_id,
#    azbin=10, win=[-50, 500], rayp=10,
#    obs_dir='obs', syn_dir='syn', syn_band_code='MX',
#    syn_suffix='.sem.sac', filtpad=100,
#    use_filter=True, freqlim=[0.009,0.011,0.05,0.1],
#    use_prefilt=True, prefilt=[0.13, 0.07],
#    use_STF=True, STF=STF)

#------
#print "\nload data\n"
#misfit.load(filename='misfit/misfit.json')
#event_id = [key for key in misfit.data['events']][0]

#------ plot misfit
#print "\nplot_misfit\n"
#for window_id in window_id_list:
#    out_file = '%s/misfit_%s.pdf' % (adj_dir, window_id)
#    misfit.plot_misfit(event_id=event_id, window_id=window_id, out_file=out_file)

#------ relocate_1d
#print "\nrelocate_1d\n"
#
#out_file = '%s/relocate_1d_fix_depth' % (adj_dir)
#reloc = misfit.relocate_1d(
#        min_cc_max=0.7,
#        min_cc_0=0.5,
#        min_SNR=10,
#        max_cc_time_shift=8.0,
#        event_id=event_id,
#        window_id_list=window_id_list,
#        fix_depth=True,
#        out_cmt_file=out_file+'.CMTSOLUTION')
#with open(out_file+'.json', 'w') as fp:
#    json.dump(reloc, fp, indent=2)
#
#out_file = '%s/relocate_1d' % (misfit_dir)
#reloc = misfit.relocate_1d(
#        min_cc_max=0.7,
#        min_cc_0=0.5,
#        min_SNR=10,
#        max_cc_time_shift=8.0,
#        event_id=event_id,
#        window_id_list=window_id_list,
#        fix_depth=False,
#        out_cmt_file=out_file+'.CMTSOLUTION')
#with open(out_file+'.json', 'w') as fp:
#    json.dump(reloc, fp, indent=2)