#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np

# read command line args
misfit_file = "misfit/misfit.pkl"
figure_dir = "misfit"

#obs_dir = str(sys.argv[1])
#syn_dir = str(sys.argv[2])
#misfit_dir = str(sys.argv[3])

#
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== load data\n"
misfit.load(filename=misfit_file)

#------
#print "\nplot misfit\n"
#
#window_id = 'F.s,S'
#out_file = "%s/%s_%s.pdf" % (misfit_dir, event_id, window_id)
#misfit.plot_misfit(event_id, 
#        window_id=window_id,
#        out_file=out_file)
#
#window_id = 'F.p,P'
#out_file = "%s/%s_%s.pdf" % (misfit_dir, event_id, window_id)
#misfit.plot_misfit(event_id, 
#        window_id=window_id,
#        out_file=out_file)

#------
print "\n====== plot seismograms\n"

plot_param = {
  'time':[-50,200], 'rayp':10., 'azbin':30, 'window_id':'F.p,P',
  'SNR':None, 'CC0':None, 'CCmax':None, 'dist':None }

misfit.plot_seismogram(
    savefig=True,
    out_dir=figure_dir,
    plot_param=plot_param)
 
#window_id = 'F.s,S'
#misfit.plot_seismograms(event_id,
#    azbin=5, win=[30, 300], rayp=16,
#    obs_dir=obs_dir,
#    syn_dir=syn_dir, syn_band_code='MX', syn_suffix='.sem.sac', 
#    out_dir=misfit_dir, savefig=True,
#    use_STF=True, 
#    use_window=True, window_id=window_id,
#    min_SNR=10, min_CCmax=0.5)