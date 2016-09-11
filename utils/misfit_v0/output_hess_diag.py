#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
model_name = str(sys.argv[2])
out_file = str(sys.argv[3])

print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== output hess diag\n")
misfit.output_hess_diag(model_name, out_file=out_file)
