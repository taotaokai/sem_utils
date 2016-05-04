#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np

# read command line args
data_dir = "DATA" 
cmt_file = "DATA/CMTSOLUTION.init"
out_file = "DATA/CMTSOLUTION.harvard"
out_file2 = "DATA/CMTSOLUTION.ecef"

#
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== setup event\n")
misfit.setup_event(cmt_file, ECEF=True)

print("\n====== write cmtsolution\n")
misfit.write_cmtsolution(out_file, ECEF=False)

print("\n====== setup event\n")
misfit.setup_event(out_file, ECEF=False)

print("\n====== write cmtsolution\n")
misfit.write_cmtsolution(out_file2, ECEF=True)