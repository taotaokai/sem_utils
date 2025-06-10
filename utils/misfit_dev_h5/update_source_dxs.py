#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("dx", type=float)
parser.add_argument("dy", type=float)
parser.add_argument("dz", type=float)

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file) 
misfit.update_source_dxs(args.dx, args.dy, args.dz)