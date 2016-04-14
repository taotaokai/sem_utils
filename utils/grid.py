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
#dm = {
#    'dt0': np.linspace(-10,0,11),
#    'dxs': np.linspace(-20,20,11),
#    }
#dm = {
#    'dt0': np.linspace(-10,0,11),
#    'dtau': np.linspace(-2,10,13),
#    }
dm = {
    'dt0': np.linspace(-10,0,11),
    'dmt': np.linspace(-50,50,11),
    }
models = [x for x in dm]
outfig = "cc_grid2_%s_%s.pdf" % (models[0], models[1])
xmax, ymax, ccmax, weight = \
  misfit.grid2_cc_perturbed_seismogram(dm=dm, plot=True, outfig=outfig, 
      plot_seism=False)

print models
print misfit.data['src_perturb']

print xmax, ymax, ccmax, weight
