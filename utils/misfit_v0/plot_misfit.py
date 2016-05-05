#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
figure_dir = str(sys.argv[2])
#misfit_file = "misfit/misfit.pkl"
#figure_dir = "misfit"

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(filename=misfit_file)

#print "\n====== plot P seismograms\n"
#plot_param = {
#  'time':[0,200], 'rayp':10., 'azbin':5, 'window_id':'F.p,P',
#  'SNR':10, 'CC0':0.5, 'CCmax':0.6, 'dist':None }
## 'SNR':-100, 'CC0':0.0, 'CCmax':0.0, 'dist':None }
#misfit.plot_seismogram(
#    savefig=True,
#    out_dir=figure_dir,
#    plot_param=plot_param)
#
#print "\n====== plot S seismograms\n"
#plot_param = {
#  'time':[0,210], 'rayp':19., 'azbin':5, 'window_id':'F.s,S',
#  'SNR':10, 'CC0':0.5, 'CCmax':0.6, 'dist':None }
## 'SNR':-100, 'CC0':0.0, 'CCmax':0.0, 'dist':None }
#misfit.plot_seismogram(
#    twopass_filt=True,
#    savefig=True,
#    out_dir=figure_dir,
#    plot_param=plot_param)

print("\n====== plot P seismograms\n")
plot_param = {
  'time':[-50,100], 'rayp':6., 'azbin':30, 'window_id':'Z.p,P',
  'SNR':10, 'CC0':0.0, 'CCmax':0.6, 'dist':None }
# 'SNR':-100, 'CC0':0.0, 'CCmax':0.0, 'dist':None }
misfit.plot_seismogram_1comp(
    savefig=True,
    out_dir=figure_dir,
    plot_param=plot_param)

plot_param = {
  'time':[-50,100], 'rayp':6., 'azbin':30, 'window_id':'R.p,P',
  'SNR':10, 'CC0':0.0, 'CCmax':0.6, 'dist':None }
# 'SNR':-100, 'CC0':0.0, 'CCmax':0.0, 'dist':None }
misfit.plot_seismogram_1comp(
    savefig=True,
    out_dir=figure_dir,
    plot_param=plot_param)

#print "\n====== plot S seismograms\n"
#plot_param = {
#  'time':[0,210], 'rayp':19., 'azbin':5, 'window_id':'F.s,S',
#  'SNR':10, 'CC0':0.5, 'CCmax':0.6, 'dist':None }
## 'SNR':-100, 'CC0':0.0, 'CCmax':0.0, 'dist':None }
#misfit.plot_seismogram(
#    twopass_filt=True,
#    savefig=True,
#    out_dir=figure_dir,
#    plot_param=plot_param)