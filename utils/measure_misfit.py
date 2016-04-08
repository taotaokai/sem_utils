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

syn_band_code = "MX"
syn_suffix = ".sem.sac"
left_pad = 100.0
right_pad = 0.0
obs_preevent = 100.0

window_list = [
    ('F','p,P', [-30,70]), 
    ('F','s,S', [-30,70]) ]
filter_param=('butter', 2, [0.01, 0.1])
taper_param=('cosine', 0.1)
weight_param={'SNR':[10, 15], 'CCmax':[0.6,0.8], 'CC0':[0.5,0.7]}

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
misfit_file = "%s/misfit.pkl" % (misfit_dir)

#
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== load data\n"
misfit.load(filename=misfit_file)

#print "\n====== setup event\n"
#misfit.setup_event(CMTSOLUTION_file, is_ECEF=True)

#print "\n====== setup station\n"
#misfit.setup_station(channel_file)

#print "\n====== read seismogram: obs, grf\n"
#misfit.read_obs_grf(
#  obs_dir=obs_dir, 
#  syn_dir=syn_dir, syn_band_code=syn_band_code, syn_suffix=syn_suffix,
#  left_pad=left_pad, right_pad=right_pad, obs_preevent=obs_preevent)

#print "\n====== setup window\n"
#misfit.setup_window(
#    window_list=window_list,
#    filter_param=filter_param,
#    taper_param=taper_param)

print "\n====== measure window\n"
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