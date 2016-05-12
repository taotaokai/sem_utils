#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from misfit import Misfit

# read command line args
misfit_file = str(sys.argv[1])

syn_band_code = "MX"
syn_suffix = ".sem.sac"

outdir_dxs = "output_dxs"
outdir_dmt = "output_dmt"

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== waveform_der_dxs\n")
misfit.waveform_der_dxs(
    syn_dir=outdir_dxs,
    syn_band_code=syn_band_code, 
    syn_suffix=syn_suffix)
#   sac_dir='du_dxs')
 
print("\n====== waveform_der_dmt\n")
misfit.waveform_der_dmt(
    syn_dir=outdir_dmt,
    syn_band_code=syn_band_code, 
    syn_suffix=syn_suffix)
#   sac_dir='du_dmt')

print("\n====== save data\n")
misfit.save(misfit_file)
