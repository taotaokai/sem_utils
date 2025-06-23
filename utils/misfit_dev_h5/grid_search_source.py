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
parser.add_argument("out_txt")  # "search.txt"
parser.add_argument("out_fig")  # "search.pdf"
parser.add_argument("-i", "--niter", default=5, type=int)  # num of iterations
parser.add_argument("-n", "--nproc", default=5, type=int)  # num of processes

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file)

print(f"======= [{datetime.now()}] grid_search_source")

misfit.grid_search_source(args.out_txt, args.out_fig, max_niter=args.niter, nproc=args.nproc)

print(f"======= [{datetime.now()}] END")