#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np
import matplotlib.pyplot as plt

# read command line args
data_dir = "DATA" 
obs_dir = "obs"
syn_dir = "output_green"
adj_dir = "adj"
misfit_dir = "misfit"
freqmin = 0.01
freqmax = 0.1
syn_band_code = "MX"
syn_suffix = ".sem.sac"
#srcfrechet_file = "output_srcfrechet/src_frechet.000001"
#outdir_dxs = "output_dxs"
#outdir_dmt = "output_dmt"

#------
print "\n====== initialize\n"
misfit = Misfit()

#------
print "\n====== load data\n"
misfit.load(filename='%s/misfit.pkl' % (misfit_dir))

print "\n====== measure_hessian_src\n"
misfit.measure_hessian_src()

print "\n====== update_source\n"
misfit.update_source()
