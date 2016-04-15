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

print "\n====== grid_cc \n"
dm = {
    't0': np.linspace(-10,0,11),
    'tau': np.linspace(-2,5,8),
    'xs': np.linspace(-10,10,21),
    }
axes = [ ['tau','xs'], ['tau','t0'], ['t0','xs'] ]
outfig = "grid_cc.pdf"

model, ccmax, weight = misfit.grid_cc_perturbed_seismogram(
    dm=dm, axes=axes, outfig=outfig, plot_seism=False)

print model, ccmax, weight
