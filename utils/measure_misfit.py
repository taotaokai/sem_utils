#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np

# read command line args
data_dir = "DATA" 
cmt_file = "DATA/CMTSOLUTION.init"
channel_file = "DATA/channel.txt"
obs_dir = "obs"
syn_dir = "output_green"
adj_dir = "adj"
misfit_file = "misfit/misfit.pkl"

sampling_rate = 1.0
before_first_arrival=100.0
after_first_arrival=200.0

syn_band_code = "MX"
syn_suffix = ".sem.sac"

window_list = [ ('F','p,P', [-30,120]) ]
filter_param=('butter', 2, [0.01, 0.08])
taper_param=('cosine', 0.1)
#weight_param={'SNR':[10, 15], 'CC0':[0.2,0.7]}
weight_param={}

#
print("\n====== initialize\n")
misfit = Misfit(sampling_rate=sampling_rate)

print("\n====== setup event\n")
misfit.setup_event(cmt_file, is_ECEF=True)

print("\n====== setup station\n")
misfit.setup_station(channel_file)

#print([x for x in misfit.data['station']]

print("\n====== read observed data \n")
misfit.read_observed_data(obs_dir,
    before_first_arrival=before_first_arrival,
    after_first_arrival=after_first_arrival)

print("\n====== read synthetic green's function\n")
misfit.read_synthetic_green(syn_dir,
    band_code=syn_band_code, suffix=syn_suffix)

print("\n====== setup window\n")
misfit.setup_window(
    window_list=window_list,
    filter_param=filter_param,
    taper_param=taper_param)

print("\n====== measure window\n")
misfit.measure_window(
    plot=False,
    weight_param=weight_param)

print("\n====== save data\n")
misfit.save(filename=misfit_file)

print("\n====== make adjoint source\n")
misfit.make_adjoint_source(
    adj_type='dchi_dg',
    out_dir=adj_dir,
    syn_band_code=syn_band_code)