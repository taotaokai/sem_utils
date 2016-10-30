#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""calculate cc for step sizes
"""
import sys
import numpy as np
import importlib.util

from misfit import Misfit

#------ read command line args
par_file = str(sys.argv[1])
misfit_file = str(sys.argv[2])
out_file = str(sys.argv[3])

#------ load parameter file
if sys.version_info < (3, ):
  raise Exception("need python3")
elif sys.version_info < (3, 5):
  spec =importlib.machinery.SourceFileLoader("misfit_par", par_file)
  par = spec.load_module()
else:
  spec = importlib.util.spec_from_file_location("misfit_par", par_file)
  par = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(par)

dt0_max = -0.057114
xs_mt_max = 0.733467

tau_1d = np.linspace(-2.0,-1.0,11)
dm_tau = {
    't0': dt0_max * np.ones(tau_1d.size),
    'xs_mt': xs_mt_max * np.ones(tau_1d.size),
    'tau': tau_1d,
    }

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== grid_cc \n")
wcc_sum, weight_sum = \
    misfit.cc_linearized_seismogram_for_source(
        dm=dm_tau, 
        plot=False)

with open(out_file, 'w') as f:
  f.write("#weight_sum = {:12.5e}\n".format(weight_sum))
  f.write("#dtau wcc_sum/weight_sum\n")
  for idx in range(len(wcc_sum)):
    f.write("{:12.5e} {:15.8e}\n".format(
      dm_tau['tau'][idx],
      wcc_sum[idx]/weight_sum))
