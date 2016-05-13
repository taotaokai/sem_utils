#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
import importlib.util
from misfit import Misfit

#------ read command line args
misfit_file = str(sys.argv[1])
par_file = str(sys.argv[2])

# output
cmt_file = str(sys.argv[3])
log_file = str(sys.argv[4])

#------ load parameter file
spec = importlib.util.spec_from_file_location("misfit_par", par_file)
par = importlib.util.module_from_spec(spec)
spec.loader.exec_module(par)

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== grid_cc \n")
#dm = {
#    't0': [-5,5],
#    'tau':[-5,5],
#    'xs': [-5,5],
#    'mt': [-5,5],
#    }

dm_opt = misfit.search1d_cc_perturbed_seismogram(
    dm_range=par.dm, 
    ngrid=par.ngrid,
    max_niter=par.max_niter,
    range_ratio=par.range_ratio,
    plot_seism=False,
    log_file=log_file,
    cmt_file=cmt_file,
    )
