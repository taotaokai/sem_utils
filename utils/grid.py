#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

import numpy as np
import matplotlib.pyplot as plt

# read command line args
misfit_file = "misfit/misfit.pkl"

#------
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== load data\n"
misfit.load(filename=misfit_file)

print "\n====== grid2 \n"
dm = {
    'dt0': np.linspace(-10,0,11),
    'dxs': np.linspace(-5,5,11),
    }
misfit.cc_perturb_grid2(dm=dm, plot=True, outfig=None)
