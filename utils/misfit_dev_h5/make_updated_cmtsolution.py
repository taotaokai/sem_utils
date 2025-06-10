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
parser.add_argument("out_cmt")         # "CMTSOLUTION.updated"
parser.add_argument("dt0", type=float)
parser.add_argument("dtau", type=float)
parser.add_argument("dxs", type=float)
parser.add_argument("dmt", type=float)

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file)
misfit.make_updated_cmtsolution(args.out_cmt,
                                dt0 = args.dt0,
                                dtau = args.dtau,
                                dxs = args.dxs,
                                dmt = args.dmt)