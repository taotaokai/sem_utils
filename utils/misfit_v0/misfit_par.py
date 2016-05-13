#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parameters
"""

#====== read_data
syn_band_code = "MX"
syn_suffix = ".sem.sac"
left_pad = 100.0
right_pad = 0.0
obs_preevent = 100.0

#====== measure_misfit
window_list = [
   ('Z','p,P', [-30,40]), 
   ('R','p,P', [-30,40]), 
   ('Z','sP', [-40,50]), 
   ('R','sP', [-40,50]), 
   ('Z','s,S', [-40,60]),
   ('R','s,S', [-40,60]),
   ('T','s,S', [-40,60]),
   ('Z','sS', [-40,80]),
   ('R','sS', [-40,80]),
   ('T','sS', [-40,80]),
   ]

filter_param=('butter', 2, [0.01, 0.08])

taper_param=('cosine', 0.1)

weight_param={
  'SNR':[10, 15],
  'CCmax':[0.6,0.8],
  'CC0':[0.5,0.7],
# 'dist':[25,26],
    }

#====== make_cmt_der
norm_dxs = 2000.0
ratio_M0 = 0.1
fix_M0 = False
zerotrace = True

#====== waveform_der
outdir_dxs = "output_dxs"
outdir_dmt = "output_dmt"

#====== search 
dm = {
    't0': [-5,5],
    'tau':[-5,5],
    'xs': [-5,5],
    'mt': [-5,5],
    }
ngrid = 4
max_niter = 5
range_ratio = 0.8