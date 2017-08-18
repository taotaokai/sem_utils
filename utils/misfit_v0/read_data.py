#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
import os
import importlib
from misfit import Misfit

# read command line args
par_file = str(sys.argv[1])
misfit_file = str(sys.argv[2])
cmt_file = str(sys.argv[3])
channel_file = str(sys.argv[4])
syn_dir = str(sys.argv[5])
obs_dir = str(sys.argv[6])

# load parameter file
if sys.version_info < (3, ):
  raise Exception("need python3")
elif sys.version_info < (3, 5):
  spec =importlib.machinery.SourceFileLoader("misfit_par", par_file)
  par = spec.load_module()
else:
  spec = importlib.util.spec_from_file_location("misfit_par", par_file)
  par = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(par)

print("\n====== initialize\n")
misfit = Misfit()

print("\n====== setup event\n")
misfit.setup_event(cmt_file, ECEF=par.cmt_is_ECEF, GPS_ELLPS=par.GPS_ELLPS)

print(misfit.data['event'])

print("\n====== setup station\n")
misfit.setup_station(channel_file, band_code=par.obs_band_code,three_channels=True, GPS_ELLPS=par.GPS_ELLPS)

print("\n====== read seismogram: obs, syn\n")
misfit.read_obs_syn(
  obs_dir=obs_dir,
  syn_dir=syn_dir, 
  syn_band_code=par.syn_band_code, 
  syn_suffix=par.syn_suffix,
  left_pad=par.left_pad, 
  right_pad=par.right_pad, 
  obs_preevent=par.obs_preevent,
  syn_is_grn=par.syn_is_grn, 
  )

print("\n====== save data\n")
misfit.save(misfit_file)