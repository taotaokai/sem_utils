#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(filename=misfit_file)

print("\n====== grid_cc \n")
dm = {
    't0': [-5,5],
    'tau':[-5,5],
    'xs': [-5,5],
#   'mt': [-5,5],
    }

dm_opt = misfit.search1d_cc_perturbed_seismogram(
    dm_range=dm, 
    ngrid=4,
    max_niter=5,
    range_ratio=0.8,
    plot_seism=False,
    log_file="search1d.log",
    cmt_file="CMTSOLUTION.search1d",
    )
