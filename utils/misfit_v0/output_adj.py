#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
adj_dir = str(sys.argv[2])

syn_band_code = "MX"

print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(filename=misfit_file)

print("\n====== output adjoint source\n")
misfit.output_adj(
    adj_type='dchi_dg',
    out_dir=adj_dir,
    syn_band_code=syn_band_code)
