#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
#import numpy as np

# read command line args
misfit_file = str(sys.argv[1])

window_list = [
   ('Z','p,P', [-30,75]), 
   ('R','p,P', [-30,75]), 
   ('Z','s,S', [-30,80]),
   ('R','s,S', [-30,80]),
   ('T','s,S', [-30,80]),
   ]
#window_list = [
#    ('F','p,P', [-30,75]),
#    ('F','s,S', [-30,80]),
#    ]

filter_param=('butter', 2, [0.01, 0.08])
taper_param=('cosine', 0.1)
weight_param={
  'SNR':[10, 15],
  'CCmax':[0.6,0.8],
  'CC0':[0.5,0.7],
# 'dist':[25,26],
    }

print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(filename=misfit_file)

print("\n====== setup window\n")
misfit.setup_window(
    window_list=window_list,
    filter_param=filter_param,
    taper_param=taper_param)

print("\n====== measure window\n")
misfit.measure_adj(
    plot=False,
    weight_param=weight_param)

print("\n====== save data\n")
misfit.save(filename=misfit_file)