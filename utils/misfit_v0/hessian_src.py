#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
#import numpy as np
#import matplotlib.pyplot as plt

# read command line args
misfit_file = "misfit/misfit.pkl"

#------
print "\n====== initialize\n"
misfit = Misfit()

#------
print "\n====== load data\n"
misfit.load(filename=misfit_file)

print "\n====== measure_hessian_src\n"
misfit.measure_hessian_src()

print "\n====== update_source\n"
misfit.update_source()
