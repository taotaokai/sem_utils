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

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file)

print(f"======= [{datetime.now()}] grid_search_source")

misfit.grid_search_structure(nproc=args.nproc)

print(f"======= [{datetime.now()}] END")