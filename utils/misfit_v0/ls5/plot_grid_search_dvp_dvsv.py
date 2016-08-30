#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from scipy import interpolate

import importlib.util

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

#------ read command line args
par_file = sys.argv[1]
fname_list = sys.argv[2]

#------ load parameter file
if sys.version_info < (3, ):
  raise Exception("need python3")
elif sys.version_info < (3, 5):
  spec =importlib.machinery.SourceFileLoader("misfit_par", par_file)
  par = spec.load_module()
else:
  spec = importlib.util.spec_from_file_location("misfit_par", par_file)
  par = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(par)

#------ read file list
with open(fname_list, "r") as f:
  # ignore blank/commeted lines
  fnames = [ l for l in (line.strip() for line in f) if l and not(l.startswith('#')) ]

evid = [ l.split('.')[0] for l in fnames ]

#------ sum up all cc list 
wcc_sum_total = np.zeros(par.vp2d.shape)
weight_sum_total = 0.0
vp_max_event = []
vsv_max_event = []
weight_sum_event = []
for fname in fnames:
  # read in data file
  with open(fname, "r") as f:
    l = f.readline().replace('\n','').split()
    weight_sum = float(l[2])
    weight_sum_event.append(weight_sum)
  with open(fname, "r") as f:
    # ignore blank/commeted lines
    lines = [ l for l in (line.strip() for line in f) if l and not(l.startswith('#')) ]
    # split each line into a list 
    lines = [ l.replace('\n','').split() for l in lines ]
  wcc = np.array([ float(x[2]) for x in lines ])
  wcc = wcc.reshape(par.vp2d.shape)
  #
  wcc_sum_total += wcc * weight_sum
  weight_sum_total += weight_sum
  #
  interp = interpolate.RectBivariateSpline(par.vsv1d, par.vp1d, wcc)
  y = np.linspace(np.min(par.vsv1d), np.max(par.vsv1d), 100)
  x = np.linspace(np.min(par.vp1d), np.max(par.vp1d), 100)
  zz = interp(y, x)
  idx_max = np.argmax(zz)
  iy, ix = np.unravel_index(idx_max, zz.shape)
  vp_max_event.append(x[ix])
  vsv_max_event.append(y[iy])

#------ plot weighted average cc
wcc2d = wcc_sum_total/weight_sum_total

#NOTE: the dimension order is reversed from np.meshgrid/plt.contour 
# x->row y->column for RectBivariateSpline
interp = interpolate.RectBivariateSpline(par.vsv1d, par.vp1d, wcc2d)
y = np.linspace(np.min(par.vsv1d), np.max(par.vsv1d), 200)
x = np.linspace(np.min(par.vp1d), np.max(par.vp1d), 200)
zz = interp(y, x)

plt.figure(figsize=(11,8.5)) # US Letter

levels = np.linspace(np.min(wcc2d), np.max(wcc2d), 100)
#cs = plt.contour(par.vp2d, par.vsv2d, wcc2d, levels=levels)
cs = plt.contour(x, y, zz, levels=levels)
plt.clabel(cs, levels[1::2], inline=1, fmt='%.3f', fontsize=10)

# mark maximum for each event
plt.plot(vp_max_event, vsv_max_event, 'ko', markersize=3, clip_on=False)
for i in range(len(vp_max_event)):
  plt.text(vp_max_event[i], vsv_max_event[i], str(evid[i]), fontsize=5,
      verticalalignment='bottom')
  plt.text(vp_max_event[i], vsv_max_event[i], str(weight_sum_event[i]), fontsize=5,
      verticalalignment='top')

# mark the maximum
idx_max = np.argmax(zz)
iy, ix = np.unravel_index(idx_max, zz.shape)
x_max = x[ix]
y_max = y[iy]
plt.plot(x_max, y_max, 'r*', markersize=10)

plt.axes().set_aspect('equal')

plt.xlabel('vp step size')
plt.ylabel('vsv step size')
plt.title('vp-vsv grid search (max vp %.3f vsv %.3f)' % (x_max, y_max))

#plt.show()
plt.savefig("grid_search_vp_vsv.pdf", format='pdf')