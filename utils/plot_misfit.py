#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import json
import numpy as np

# read command line args
obs_dir = str(sys.argv[1])
syn_dir = str(sys.argv[2])
misfit_dir = str(sys.argv[3])

#------
print "\ninitialize\n"
misfit = Misfit()

#------
print "\nload data\n"
misfit.load(filename='%s/misfit.json' % (misfit_dir))

event_id = [ key for key in misfit.data['events'] ][0]

#------
print "\nplot misfit\n"

window_id = 'F.s,S'
out_file = "%s/%s_%s.pdf" % (misfit_dir, event_id, window_id)
misfit.plot_misfit(event_id, 
        window_id=window_id,
        out_file=out_file)

window_id = 'F.p,P'
out_file = "%s/%s_%s.pdf" % (misfit_dir, event_id, window_id)
misfit.plot_misfit(event_id, 
        window_id=window_id,
        out_file=out_file)

#------
print "\nplot seismograms\n"

window_id = 'F.p,P'
misfit.plot_seismograms(event_id,
    azbin=5, win=[0, 150], rayp=10,
    obs_dir=obs_dir, 
    syn_dir=syn_dir, syn_band_code='MX', syn_suffix='.sem.sac', 
    out_dir=misfit_dir, savefig=True,
    use_STF=True, 
    use_window=True, window_id=window_id,
    min_SNR=10, min_CCmax=0.5)
 
window_id = 'F.s,S'
misfit.plot_seismograms(event_id,
    azbin=5, win=[30, 300], rayp=16,
    obs_dir=obs_dir,
    syn_dir=syn_dir, syn_band_code='MX', syn_suffix='.sem.sac', 
    out_dir=misfit_dir, savefig=True,
    use_STF=True, 
    use_window=True, window_id=window_id,
    min_SNR=10, min_CCmax=0.5)