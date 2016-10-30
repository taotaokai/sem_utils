#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
model_name = str(sys.argv[2])
adj_dir = str(sys.argv[3])

syn_band_code = "MX"

print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== output adjoint source\n")
misfit.output_adj_hess_model_product(
    model_name=model_name,
    out_dir=adj_dir,
    syn_band_code=syn_band_code)
