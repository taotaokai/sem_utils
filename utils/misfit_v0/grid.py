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
    't0': np.linspace(-1,1,4),
#   'tau': np.linspace(-2,2,5),
    'xs': np.linspace(1.2, 2.2, 4),
    'mt': np.linspace(1.39,1.39,1),
    }
#axes = [ ['t0', 'tau'], ['t0','xs'], ['t0', 'mt'], 
#         ['tau','xs'],  ['tau', 'mt'], ['xs', 'mt'],
#        ['t0',], ['tau',], ['xs',], ['mt',]
#         ]

#axes = [ ['mt','t0'] ]
#axes = [ ['tau','t0'] ]
axes = [ ['xs', 't0'] ]
#axes = [ ['xs', 'mt'] ]
#axes = [ [x for x in dm] ]
#axes = [ ['xs'], ['t0'], ['xs','t0'] ]

outfig = "grid_cc_%s.pdf" % ("_".join([x for x in dm]))
#outfig = "grid_cc.pdf"

model, ccmax, weight = misfit.grid_cc_perturbed_seismogram(
    dm=dm, axes=axes, outfig=outfig, plot_seism=False,
    cmt_file="CMTSOLUTION.grid_cc"
    )

print model, ccmax, weight
