#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parameters
"""
import numpy as np

#====== read_data
syn_band_code = "MX"
syn_suffix = ".sem.sac"
left_pad = 100.0
right_pad = 0.0
obs_preevent = 100.0
syn_is_grn = False
cmt_is_ECEF = True

#====== measure_misfit
# misfit window specfications
def make_window_list_P_wave(evdp_km):
  flo = 0.01
  fhi = 0.05
  if evdp_km <= 150:
    window_list_P_wave = [
        {'phase':'p,P,pP,sP', 'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'p,P,pP,sP', 'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'PP',        'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[25,180]},
        {'phase':'PP',        'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[25,180]},
        ]
  elif evdp_km > 150 and evdp_km <=400:
    window_list_P_wave = [
        {'phase':'p,P',      'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'p,P',      'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'pP,sP,PP', 'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'pP,sP,PP', 'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        ]
  else: # evdp_km > 400:
    window_list_P_wave = [
        {'phase':'p,P',   'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'p,P',   'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'pP,PP', 'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'pP,PP', 'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sP',    'component':'Z', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sP',    'component':'R', 'time':[-30,50], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        ]
  return window_list_P_wave

def make_window_list_S_wave(evdp_km):
  flo = 0.01
  fhi = 0.05
  if evdp_km <= 150:
    window_list_S_wave = [
        {'phase':'s,S,sS', 'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'s,S,sS', 'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'s,S,sS', 'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'SS',     'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[28,180]},
        {'phase':'SS',     'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[28,180]},
        {'phase':'SS',     'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[28,180]},
        ]
  elif evdp_km > 150 and evdp_km <=300:
    window_list_S_wave = [
        {'phase':'s,S', 'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'s,S', 'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'s,S', 'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'SS',  'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[30,180]},
        {'phase':'SS',  'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[30,180]},
        {'phase':'SS',  'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[30,180]},
        {'phase':'sS',  'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sS',  'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sS',  'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        ] 
  else:
    window_list_S_wave = [
        {'phase':'s,S',   'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'s,S',   'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'s,S',   'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sS,SS', 'component':'Z', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sS,SS', 'component':'R', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        {'phase':'sS,SS', 'component':'T', 'time':[-40,60], 'filter':[flo, fhi, 2], 'dist':[0,180]},
        ]
  return window_list_S_wave

def make_window_list_surface_wave(evdp_km):
  # for 50s period, the phase velocity sensitivity kernel can only reach to about 200/150 km depth for fundamental Rayleigh/Love waves, respectively  
  # slowness: 22 ~ 44 s/deg = 5 ~ 2.5 km/s
  flo = 0.01
  fhi = 0.025
  if evdp_km <= 150:
    window_list_surface_wave = [
       {'phase':'surface', 'component':'Z', 'time':[-50,50], 'slowness':[22,44], 'filter':[flo, fhi, 2]},
       {'phase':'surface', 'component':'R', 'time':[-50,50], 'slowness':[22,44], 'filter':[flo, fhi, 2]},
       {'phase':'surface', 'component':'T', 'time':[-50,50], 'slowness':[22,44], 'filter':[flo, fhi, 2]},
       ]
  else:
    window_list_surface_wave = []

  return window_list_surface_wave

# misfit window weight scheme
weight_param={
  'cc_tshift':[-10,-8,8,10],
  'SNR':[10, 15],
  'CCmax':[0.6,0.8],
  'CC0':[0.5,0.7],
# 'dist':[25,26],
    }

#====== plot configurations

# plot_seismogram_1comp
plot_azbin_size = 10
plot_begin_time = -50
plot_end_time = 50
plot_clip_ratio = 2.0

#====== source inversion
# make_cmt_der
norm_dxs = 2000.0
ratio_M0 = 0.1
fix_M0 = False
zerotrace = True

# waveform_der
outdir_dxs = "output_dxs"
outdir_dmt = "output_dmt"

# search 
dm = {
    't0': [-5,5],
    'tau':[-5,5],
    'xs': [-5,5],
    'mt': [-5,5],
    }
ngrid = 4
max_niter = 5
range_ratio = 0.8

#====== structure inversion
# search range
dmodel_step_size = np.arange(0, 10, 0.1)
