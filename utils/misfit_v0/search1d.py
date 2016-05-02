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

dist_min=0.0

#------
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== load data\n"
misfit.load(filename=misfit_file)

print "\n====== grid_cc \n"
dm = {
    't0': [-2,2],
    'tau':[-2,2],
    'xs': [-5,5],
    'mt': [-5,5],
    }

misfit.search1d_cc_perturbed_seismogram(
    dm_range=dm, 
    ngrid=4,
    range_ratio=0.85,
    dist_min=dist_min,
    plot_seism=True,
    )
