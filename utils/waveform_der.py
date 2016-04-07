#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import numpy as np
import matplotlib.pyplot as plt

# read command line args
data_dir = "DATA" 
obs_dir = "obs"
syn_dir = "output_green"
adj_dir = "adj"
misfit_dir = "misfit"
freqmin = 0.01
freqmax = 0.1
syn_band_code = "MX"
syn_suffix = ".sem.sac"
srcfrechet_file = "output_adj/src_frechet.000001"

#------
print "\n====== initialize\n"
misfit = Misfit()

#------
print "\n====== load data\n"
misfit.load(filename='%s/misfit.pkl' % (misfit_dir))

#------
print "\n====== read srcfrechet\n"
misfit.read_srcfrechet(filename=srcfrechet_file, update=True)

#------
print "\n====== make_cmt_dxs/dmt\n"
out_file='%s/CMTSOLUTION.dxs'%(data_dir)
misfit.make_cmt_dxs(out_file=out_file, norm=2000.0)

out_file='%s/CMTSOLUTION.dmt'%(data_dir)
misfit.make_cmt_dmt(out_file=out_file, fix_M0=True)

#------
print "\n====== save data\n"
misfit.save(filename='%s/misfit.pkl' % (misfit_dir))

#mt = misfit.data['event']['mt']
#dmt = misfit.data['src_perturb']['dmt']
#
#print mt, dmt
#print np.sum(dmt*mt)/np.sum(dmt**2)**0.5/np.sum(mt**2)**0.5

#------
#print "\n====== waveform_der_stf\n"
#misfit.waveform_der_stf()
#
#print "\n====== waveform_der_dxs\n"
#misfit.waveform_der_dxs(
#    syn_dir=outdir_dxs,
#    syn_band_code=syn_band_code, 
#    syn_suffix=syn_suffix)
#    #sac_dir='du_dxs')
#
#print "\n====== waveform_der_dmt\n"
#misfit.waveform_der_dmt(
#    syn_dir=outdir_dmt,
#    syn_band_code=syn_band_code, 
#    syn_suffix=syn_suffix)
#    #sac_dir='du_dmt')
