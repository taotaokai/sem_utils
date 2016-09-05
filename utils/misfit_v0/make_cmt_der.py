#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])
srcfrechet_file = str(sys.argv[2])
out_dir = str(sys.argv[3])

#misfit_file = "misfit/misfit.pkl"
#srcfrechet_file = "output_srcfrechet/src_frechet.000001"
#out_dir = "DATA"

norm_dxs = 2000.0
ratio_M0 = 0.1

#======
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== read srcfrechet\n")
misfit.read_srcfrechet(filename=srcfrechet_file)

print("\n====== make_cmt_dxs/dmt\n")
out_file='%s/CMTSOLUTION.dxs' % (out_dir)
misfit.make_cmt_dxs(out_file=out_file, norm=norm_dxs)

out_file='%s/CMTSOLUTION.dmt'%(out_dir)
misfit.make_cmt_dmt(out_file=out_file, fix_M0=False, zerotrace=True, ratio_M0=ratio_M0)
 
print("\n====== save data\n")
misfit.save(misfit_file)
