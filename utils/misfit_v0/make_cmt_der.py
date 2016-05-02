#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

# read command line args
data_dir = "DATA"
misfit_file = "misfit/misfit.pkl"
srcfrechet_file = "output_srcfrechet/src_frechet.000001"
norm_dxs = 2000.0
ratio_M0 = 0.1

#------
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== load data\n"
misfit.load(filename=misfit_file)

print "\n====== read srcfrechet\n"
misfit.read_srcfrechet(filename=srcfrechet_file, update=True)

print "\n====== make_cmt_dxs/dmt\n"
out_file='%s/CMTSOLUTION.dxs'%(data_dir)
misfit.make_cmt_dxs(out_file=out_file, norm=norm_dxs)

out_file='%s/CMTSOLUTION.dmt'%(data_dir)
misfit.make_cmt_dmt(out_file=out_file, fix_M0=True,
    zerotrace=True, ratio_M0=ratio_M0)

print "\n====== save data\n"
misfit.save(filename=misfit_file)
