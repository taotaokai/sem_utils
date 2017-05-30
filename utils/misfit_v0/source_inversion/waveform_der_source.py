#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import importlib.util
from misfit import Misfit

# read command line args
par_file = str(sys.argv[1])
misfit_file = str(sys.argv[2])
syn_dir = str(sys.argv[3])
model_name = str(sys.argv[4])

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

print("\n====== waveform_der_source\n")
misfit.waveform_der_source(
    syn_dir=syn_dir,
    syn_band_code=par.syn_band_code, 
    syn_suffix=par.syn_suffix,
    model_name=model_name,
    )
 
print("\n====== save data\n")
misfit.save(misfit_file)
