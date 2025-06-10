#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("srcfrechet_file" )  # "srcfrechet"
parser.add_argument("out_diff_cmtsolution") # diff_CMTSOLUTION
parser.add_argument("out_cmtsolution_dxs") # CMTSOLUTION_dxs
parser.add_argument("out_cmtsolution_dmt") # CMTSOLUTION_dmt

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file)
misfit.read_srcfrechet(args.srcfrechet_file)
misfit.make_perturbed_cmtsolution(args.out_diff_cmtsolution,
                                  args.out_cmtsolution_dxs,
                                  args.out_cmtsolution_dmt)