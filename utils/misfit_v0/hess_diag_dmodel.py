#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
model_name = str(sys.argv[2])

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== hess_diag_dmodel\n")
misfit.hess_diag_dmodel(model_name)
 
print("\n====== save data\n")
misfit.save(misfit_file)
