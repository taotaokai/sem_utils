#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("config_file", help="config file, e.g. misfit.yaml")

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file) 
misfit.read_config_file(args.config_file)