#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot synthetic and observed seismograms
"""
import sys
import importlib.util
from misfit import Misfit

# read command line args
par_file = str(sys.argv[1])
misfit_file = str(sys.argv[2])
figure_dir = str(sys.argv[3])

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

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== plot seismograms\n")

evdp_km = misfit.data['event']['depth']
print("evdp_km = %f" % (evdp_km))

window_list_P_wave = par.make_window_list_P_wave(evdp_km)
print(window_list_P_wave)

for win in window_list_P_wave:
  window_id = "%s_%s" % (win['phase'], win['component'])
  print("\n------ %s\n" % (window_id))
  misfit.plot_seismogram_1comp(
      savefig=True,
      out_dir=figure_dir,
      window_id=window_id,
      azbin=10,
      begin_time=-50,
      end_time=100,
      clip_ratio=2.0,
      min_CC0=min(par.weight_param['CC0']),
      min_CCmax=min(par.weight_param['CCmax']),
      min_SNR=min(par.weight_param['SNR']),
      )

window_list_S_wave = par.make_window_list_S_wave(evdp_km)
print(window_list_S_wave)

for win in window_list_S_wave:
  window_id = "%s_%s" % (win['phase'], win['component'])
  print("\n------ %s\n" % (window_id))
  misfit.plot_seismogram_1comp(
      savefig=True,
      out_dir=figure_dir,
      window_id=window_id,
      azbin=10,
      begin_time=-50,
      end_time=100,
      clip_ratio=2.0,
      min_CC0=min(par.weight_param['CC0']),
      min_CCmax=min(par.weight_param['CCmax']),
      min_SNR=min(par.weight_param['SNR']),
      )

window_list_surface_wave = par.make_window_list_surface_wave(evdp_km)
print(window_list_surface_wave)

for win in window_list_surface_wave:
  window_id = "%s_%s" % (win['phase'], win['component'])
  print("\n------ %s\n" % (window_id))
  misfit.plot_seismogram_1comp(
      savefig=True,
      out_dir=figure_dir,
      window_id=window_id,
      azbin=10,
      begin_time=-50,
      end_time=200,
      clip_ratio=2.0,
      min_CC0=min(par.weight_param['CC0']),
      min_CCmax=min(par.weight_param['CCmax']),
      min_SNR=min(par.weight_param['SNR']),
      )