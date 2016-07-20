#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot histograms of dt,cc from misfit.txt

misfit.txt: 
#station window weight CC0 CCmax dt_cc SNR AR0 ARmax
...
"""
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# read command line args
misfit_file = str(sys.argv[1])
out_file = str(sys.argv[2])

# read in misfit_file
with open(misfit_file, 'r') as f:
  lines = [x.replace('\n','').split() for x in f.readlines() if not(x.startswith('#'))]

#----- create figure
fig = plt.figure(figsize=(8.5, 11)) # US Letter

# a matrix of sub plot 
nrow = 4
ncol = 2
subplot_size = np.array([1.0/ncol, 1.0/nrow])

# axis position relative to the subplot region
ax_origin_subplot = np.array([0.2, 0.2])
# size: [width, height]
ax_size_subplot = np.array([0.7, 0.7])

# hist parameters
nbins = 50
max_dt = 10 

#----- surface RZ: CC
nrow = 1
ncol = 1
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of CC
cc = [ float(x[3]) for x in lines if (x[1]=='surface_R' or x[1]=='surface_Z') ]
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of window')
ax.set_xlim([0.5, 1])
title_str = "surface_RZ %.3f$\pm$%.3f" % (np.mean(cc), np.std(cc))
ax.set_title(title_str)

#----- surface RZ: dt
nrow = 1
ncol = 2
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of dt_cc
dt = np.array([ -1.0*float(x[5]) for x in lines if (x[1]=='surface_R' or x[1]=='surface_Z') ])
ax.hist(dt[abs(dt)<=max_dt], nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of window')
ax.set_xlim([-max_dt, max_dt])
title_str = "surface_RZ %.3f$\pm$%.3f" % (np.mean(dt), np.std(dt))
ax.set_title(title_str)

#----- surface T: CC
nrow = 2
ncol = 1
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of CC
cc = [ float(x[3]) for x in lines if (x[1]=='surface_T') ]
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of window')
ax.set_xlim([0.5, 1])
title_str = "surface_T %.3f$\pm$%.3f" % (np.mean(cc), np.std(cc))
ax.set_title(title_str)

#----- surface T: dt
nrow = 2
ncol = 2
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of dt_cc
dt = np.array([ -1.0*float(x[5]) for x in lines if (x[1]=='surface_T') ])
ax.hist(dt[abs(dt)<=max_dt], nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of window')
ax.set_xlim([-max_dt, max_dt])
title_str = "surface_T %.3f$\pm$%.3f" % (np.mean(dt), np.std(dt))
ax.set_title(title_str)

#----- body waves: CC
nrow = 3
ncol = 1
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of CC
cc = [ float(x[3]) for x in lines if ('surface' not in x[1]) ]
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of window')
ax.set_xlim([0.5, 1])
title_str = "body waves %.3f$\pm$%.3f" % (np.mean(cc), np.std(cc))
ax.set_title(title_str)

#----- body waves: dt
nrow = 3
ncol = 2
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of dt_cc
dt = np.array([ -1.0*float(x[5]) for x in lines if ('surface' not in x[1]) ])
ax.hist(dt[abs(dt)<=max_dt], nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of window')
ax.set_xlim([-max_dt, max_dt])
title_str = "body waves %.3f$\pm$%.3f" % (np.mean(dt), np.std(dt))
ax.set_title(title_str)

#----- all: CC
nrow = 4
ncol = 1
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of CC
cc = [ float(x[3]) for x in lines ]
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of window')
ax.set_xlim([0.5, 1])
title_str = "All %.3f$\pm$%.3f" % (np.mean(cc), np.std(cc))
ax.set_title(title_str)

#----- all: dt
nrow = 4
ncol = 2
# create axis
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
# get arrays of dt_cc
dt = np.array([ -1.0*float(x[5]) for x in lines ])
ax.hist(dt[abs(dt)<=max_dt], nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of window')
ax.set_xlim([-max_dt, max_dt])
title_str = "All %.3f$\pm$%.3f" % (np.mean(dt), np.std(dt))
ax.set_title(title_str)

#------ save figure
plt.savefig(out_file, format='pdf')