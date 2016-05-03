#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np

# read command line args
#data_dir = "DATA" 
#cmt_file = "DATA/CMTSOLUTION.init"
#channel_file = "DATA/channel.txt"
misfit_file = "misfit/misfit.pkl"

#
print("\n====== initialize\n")
misfit = Misfit()

#print("\n====== setup event\n")
#misfit.setup_event(cmt_file, is_ECEF=True)
#
#print("\n====== setup station\n")
#misfit.setup_station(channel_file)

print("\n====== load data\n")
misfit.load(filename=misfit_file)

print("\n====== plot radiation pattern\n")
misfit.radiation_pattern()