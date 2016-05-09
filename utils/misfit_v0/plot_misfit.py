#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot synthetic and observed seismograms
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
figure_dir = str(sys.argv[2])

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(filename=misfit_file)

print("\n====== plot P seismograms\n")
for comp in ['Z', 'R']:
  window_id = "%s.p,P" % (comp)
  plot_param = {
    'time':[-50,250], 'rayp':10., 'azbin':5, 'window_id':window_id,
    'SNR':10, 'CC0':0.5, 'CCmax':0.6, 'dist':None }
  misfit.plot_seismogram_1comp(
      savefig=True,
      out_dir=figure_dir,
      plot_param=plot_param)

print("\n====== plot S seismograms\n")
for comp in ['Z', 'R', 'T']:
  window_id = "%s.s,S" % (comp)
  plot_param = {
    'time':[-50,450], 'rayp':18., 'azbin':5, 'window_id':window_id,
    'SNR':10, 'CC0':0.5, 'CCmax':0.6, 'dist':None }
  misfit.plot_seismogram_1comp(
      savefig=True,
      out_dir=figure_dir,
      plot_param=plot_param)