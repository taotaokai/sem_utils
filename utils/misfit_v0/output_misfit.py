#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
out_file = str(sys.argv[2])

print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== output misfit measurments\n")
misfit.output_misfit(out_file)
