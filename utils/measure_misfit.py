#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np

# read command line args

data_dir = "DATA" 
obs_dir = "obs"
syn_dir = "output_green"
adj_dir = "adj"
misfit_dir = "misfit"
freqmin = 0.01
freqmax = 0.1
syn_band_code = "MX"

#data_dir = str(sys.argv[1])
#obs_dir = str(sys.argv[2])
#syn_dir = str(sys.argv[3])
#adj_dir = str(sys.argv[4])
#misfit_dir = str(sys.argv[5])
#freqmin = float(sys.argv[6])
#freqmax = float(sys.argv[7])

#------ input files
CMTSOLUTION_file = '%s/CMTSOLUTION.init' % (data_dir)
channel_file = '%s/channel.txt' % (data_dir)

#
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== setup event\n"
misfit.setup_event(CMTSOLUTION_file, is_ECEF=True)

print "\n====== setup station\n"
misfit.setup_station(channel_file)

print "\n====== read seismogram: obs, grf\n"
misfit.read_obs_grf(
  obs_dir=obs_dir, 
  syn_dir=syn_dir, syn_band_code=syn_band_code, syn_suffix=".sem.sac",
  left_pad=100, right_pad=0)

print "\n====== setup window\n"
window_list = [ 
    ('F','p,P',[-30,50]), 
]
print "window_list= ", window_list
filter_param=('butter', 2, [freqmin, freqmax])
print "filter_param= ", filter_param
misfit.setup_window(window_list=window_list,
        filter_param=filter_param)

print "\n====== measure window\n"
#weight_param={'SNR':[10, 15], 'CCmax':[0.6,0.8], 'CC0':[0.5,0.7]}
weight_param={}
print "weight_param= ", weight_param
# 
misfit.measure_adj(
    plot=False,
    weight_param=weight_param)

print "\n====== save data\n"
misfit.save(filename='%s/misfit.pkl' % (misfit_dir))

print "\n====== output adjoint source\n"
misfit.output_adj(
    adj_type='dchi_dg',
    out_dir=adj_dir,
    syn_band_code=syn_band_code)