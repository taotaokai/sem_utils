#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""calculate cc for step sizes
"""
import sys
import numpy as np
import argparse
from datetime import datetime
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("-n", "--nproc", default=5, type=int)  # num of processes
parser.add_argument("--dm_tags", nargs="+", default=['dm'])  # dmodel tags
parser.add_argument("--dm_steps", nargs="+", default=['0,4,20'])  # dmodel steps

args = parser.parse_args()

assert(len(args.dm_tags) == len(args.dm_steps))

dm = {}
for tag, val in zip(args.dm_tags, args.dm_steps):
    vals = [float(x) for x in val.split(',')]
    assert(len(vals) == 3)
    steps = np.linspace(vals[0], vals[1], vals[2])
    dm[tag] = steps

misfit = Misfit(args.misfit_h5file)

misfit.grid_search_structure(dm=dm, nproc=args.nproc)