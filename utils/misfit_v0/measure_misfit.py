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

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== delete all misfit windows\n")
misfit.delete_window()

print("\n====== add P wave windows\n")
evdp_km = misfit.data['event']['depth']
print("evdp_km = %f" % (evdp_km))

window_list_P_wave = par.make_window_list_P_wave(evdp_km)
print(window_list_P_wave)

for win in window_list_P_wave:
  misfit.add_window_body_wave(
      component=win['component'],
      phase=win['phase'],
      begin_time=win['time'][0],
      end_time=win['time'][1],
      min_frequency=win['filter'][0],
      max_frequency=win['filter'][1],
      filter_order=win['filter'][2],
      min_dist=win['dist'][0],
      max_dist=win['dist'][1],
      pre_weight=win['pre_weight'],
      )

print("\n====== add S wave windows\n")
window_list_S_wave = par.make_window_list_S_wave(evdp_km)
print(window_list_S_wave)

for win in window_list_S_wave:
  misfit.add_window_body_wave(
      component=win['component'],
      phase=win['phase'],
      begin_time=win['time'][0],
      end_time=win['time'][1],
      min_frequency=win['filter'][0],
      max_frequency=win['filter'][1],
      filter_order=win['filter'][2],
      min_dist=win['dist'][0],
      max_dist=win['dist'][1],
      pre_weight=win['pre_weight'],
      )

print("\n====== add surface wave windows\n")
window_list_surface_wave = par.make_window_list_surface_wave(evdp_km)
print(window_list_surface_wave)

for win in window_list_surface_wave:
  misfit.add_window_surface_wave(
      phase=win['phase'],
      component=win['component'],
      begin_time=win['time'][0],
      end_time=win['time'][1],
      min_slowness=win['slowness'][0],
      max_slowness=win['slowness'][1],
      min_frequency=win['filter'][0],
      max_frequency=win['filter'][1],
      filter_order=win['filter'][2],
      pre_weight=win['pre_weight'],
      )

print("\n====== measure window\n")
print("misfit_type = %s" % (par.misfit_type))
misfit.measure_adj(
    plot=False,
    misfit_type=par.misfit_type,
    weight_param=par.weight_param)

print("\n====== save data\n")
misfit.save(misfit_file)