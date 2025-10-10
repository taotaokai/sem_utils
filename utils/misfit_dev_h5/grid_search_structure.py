#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""calculate cc for step sizes
"""
import sys
import argparse
from datetime import datetime
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("-n", "--nproc", default=5, type=int)  # num of processes
parser.add_argument("--par_file", default=None)  # "misfit.h5"

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file)

if args.par_file is not None:
    misfit.read_config_file(args.par_file)

misfit.grid_search_structure(nproc=args.nproc)