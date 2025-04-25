#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""calculate cc for step sizes
"""
import sys
import warnings
import importlib.util
import datetime
import numpy as np 
from scipy import interpolate

import matplotlib
matplotlib.use("pdf")
from matplotlib.backends.backend_pdf import PdfPages 
import matplotlib.pyplot as plt

from misfit import Misfit

#------ read command line args
par_file = str(sys.argv[1])
misfit_file = str(sys.argv[2])
out_file = str(sys.argv[3])
out_figure = str(sys.argv[4])

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

#------
print("\n====== initialize\n")
misfit = Misfit()

print("\n====== load data\n")
misfit.load(misfit_file)

print("\n====== grid_cc \n")
range_shrink_ratio = 0.6
dtau_range = 2
dt0_range = 5
xs_mt_step_range = 5

dtau_opt = 0.0
dt0_opt = 0.0
xs_mt_step_opt = 0.0

niter = 0
tau0 = misfit.data['event']['tau']
#if tau0 <= 0.0:
#  warnings.warn("replace zero tau0 to 1.0!")
#  tau0 = 1.0
f = open(out_file, 'w')
with PdfPages(out_figure) as pdf:
  while niter < 10:
    f.write("====== iter = {:02d}\n".format(niter))

    # define search range
    dtau_min = dtau_opt - dtau_range
    if (tau0 + dtau_min) < 0.1: # in case of negative tau
      dtau_min = 0.1 - tau0
    dtau_max = dtau_opt + dtau_range
    dt0_min = dt0_opt - dt0_range
    dt0_max = dt0_opt + dt0_range
    xs_mt_step_min = xs_mt_step_opt - xs_mt_step_range
    xs_mt_step_max = xs_mt_step_opt + xs_mt_step_range

    f.write("search range for dtau: {:f} {:f}\n".format(dtau_min, dtau_max))
    f.write("search range for dt0:  {:f} {:f}\n".format(dt0_min, dt0_max))
    f.write("search range for xs_mt_step: {:f} {:f}\n".format(xs_mt_step_min, xs_mt_step_max))
  
    # 2D grid search for t0 and xs_mt
    dt0_1d = np.linspace(dt0_min, dt0_max, 5)
    xs_mt_step_1d = np.linspace(xs_mt_step_min, xs_mt_step_max, 5)
    dt0_2d, xs_mt_step_2d = np.meshgrid(dt0_1d, xs_mt_step_1d)
    dm = {
        'tau': dtau_opt*np.ones(dt0_2d.size),
        't0': dt0_2d.reshape(dt0_2d.size),
        'xs_mt': xs_mt_step_2d.reshape(xs_mt_step_2d.size),
        }
    wcc_sum, weight_sum = \
        misfit.cc_linearized_seismogram_for_source(
            dm=dm,
            plot=False)
    wcc = wcc_sum/weight_sum
    ind = np.argmax(wcc)
    wcc_max = wcc[ind]
    #iy, ix = np.unravel_index(ind, dt0_2d.shape)
    dt0_opt = dm['t0'][ind]
    xs_mt_step_opt = dm['xs_mt'][ind]

    #interp = interpolate.RectBivariateSpline(xs_mt_step_1d, dt0_1d, wcc)
    #y = np.linspace(np.min(xs_mt_step_1d), np.max(xs_mt_step_1d), 500)
    #x = np.linspace(np.min(dt0_1d), np.max(dt0_1d), 500)
    #zz = interp(y, x)
    #iy, ix = np.unravel_index(np.argmax(zz), zz.shape)
    #dt0_opt = x[ix]
    #xs_mt_step_opt = y[iy]

    plt.figure(figsize=(6,12)) # US Letter
    plt.subplot(2,1,1)
    #levels = np.linspace(np.min(zz), np.max(zz), 100)
    levels = np.linspace(np.min(wcc), np.max(wcc), 100)
    #cs = plt.contour(x, y, zz, levels=levels)
    cs = plt.contour(dt0_1d, xs_mt_step_1d, wcc.reshape(dt0_2d.shape), levels=levels)
    plt.clabel(cs, levels, inline=1, fmt='%.3f', fontsize=10)
    plt.plot(dm['t0'], dm['xs_mt'], 'ko', clip_on=False)
    for i in range(dt0_2d.size):
      plt.text(dm['t0'][i], dm['xs_mt'][i], '%.3f'%(wcc[i]))
    plt.plot(dt0_opt, xs_mt_step_opt, 'r+', markersize=8, clip_on=False)
    #plt.text(dt0_opt, xs_mt_step_opt, "%.3f/%.3f" % (dt0_opt, xs_mt_step_opt))
    plt.xlabel('dt0 (second)')
    plt.ylabel('xs_mt_step')
    plt.title("iter%02d dt0_opt=%f,xs_mt_step_opt=%f\n" % (niter, dt0_opt, xs_mt_step_opt))
  
    # 1D grid search for tau
    dtau_1d = np.linspace(dtau_min, dtau_max, 11)
    dm = {
        'tau': dtau_1d,
        't0': dt0_opt*np.ones(dtau_1d.size),
        'xs_mt': xs_mt_step_opt*np.ones(dtau_1d.size),
        }
    wcc_sum, weight_sum = \
        misfit.cc_linearized_seismogram_for_source(
            dm=dm,
            plot=False)
    wcc = wcc_sum/weight_sum
    #interp = interpolate.interp1d(dtau_1d, wcc, kind='cubic')
    #x = np.linspace(dtau_min, dtau_max, 500)
    #z = interp(x)
    #idx_max = np.argmax(z)
    #dtau_opt = x[idx_max]
    ind = np.argmax(wcc)
    dtau_opt = dm['tau'][ind]

    f.write("wcc_max = {:f}\n".format(wcc[ind]))
    f.write("dt0_opt = {:f}\n".format(dt0_opt))
    f.write("xs_mt_step_opt = {:f}\n".format(xs_mt_step_opt))
    f.write("dtau_opt = {:f}\n".format(dtau_opt))
  
    plt.subplot(2,1,2)
    plt.plot(dtau_1d, wcc, 'ko')
    #plt.plot(x, z, 'k')
    #plt.plot(x[idx_max], z[idx_max], 'ro', markersize=5)
    plt.plot(dtau_opt, wcc[ind], 'r+', markersize=8)
    for i in range(dtau_1d.size):
      plt.text(dm['tau'][i], wcc[i], '%.3f'%(wcc[i]))
    plt.xlabel('dtau (second)')
    plt.ylabel('wcc')
    plt.title("dau_opt=%f" % (dtau_opt))

    pdf.savefig()
    plt.close()
  
    # shrink search range
    dtau_range = range_shrink_ratio*dtau_range
    dt0_range = range_shrink_ratio*dt0_range
    xs_mt_step_range = range_shrink_ratio*xs_mt_step_range
    niter += 1

#====== close file
f.close()
