#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np

# read command line args
cmt_file = str(sys.argv[1])
out_file = str(sys.argv[2])

#
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== setup event\n")
misfit.setup_event(cmt_file, ECEF=False)

print("\n====== write cmtsolution\n")
misfit.write_cmtsolution(out_file, ECEF=True)

#print("\n====== setup event\n")
#misfit.setup_event(out_file, ECEF=False)
 
#print("\n====== write cmtsolution\n")
#misfit.write_cmtsolution(out_file2, ECEF=True)