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
matplotlib.use("pdf")
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

#----- surface RZ

# get arrays of dt(syn-obs), cc 
data_type = "surface_RZ"
data = [ [ float(x[2]), float(x[3]), -1.0*float(x[5]) ] for x in lines if (x[1]=='surface_R' or x[1]=='surface_Z') ]

weight = np.array([ x[0] for x in data ])
cc = np.array([ x[1] for x in data ])
dt = np.array([ x[2] for x in data ])

idx = np.abs(dt) <= max_dt
dt = dt[idx]
weight = weight[idx]
cc = cc[idx]

# plot CC
nrow = 1
ncol = 1
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of windows')
ax.set_xlim([0.5, 1])
title_str = "%s %.4f(%.2f)" % (data_type, np.sum(cc*weight)/np.sum(weight), np.sum(weight))
ax.set_title(title_str)

# plot dt
nrow = 1
ncol = 2
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(dt, nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of windows')
ax.set_xlim([-max_dt, max_dt])
title_str = "%s %.3f$\pm$%.3f" % (data_type, np.mean(dt), np.std(dt))
ax.set_title(title_str)


#----- surface T

# get arrays of dt(syn-obs), cc 
data_type = "surface_T"
data = [ [ float(x[2]), float(x[3]), -1.0*float(x[5]) ] for x in lines if x[1]=='surface_T' ]

weight = np.array([ x[0] for x in data ])
cc = np.array([ x[1] for x in data ])
dt = np.array([ x[2] for x in data ])

idx = np.abs(dt) <= max_dt
dt = dt[idx]
weight = weight[idx]
cc = cc[idx]

# plot CC
nrow = 2
ncol = 1
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of windows')
ax.set_xlim([0.5, 1])
title_str = "%s %.4f(%.2f)" % (data_type, np.sum(cc*weight)/np.sum(weight), np.sum(weight))
ax.set_title(title_str)

# plot dt
nrow = 2
ncol = 2
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(dt, nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of windows')
ax.set_xlim([-max_dt, max_dt])
title_str = "%s %.3f$\pm$%.3f" % (data_type, np.mean(dt), np.std(dt))
ax.set_title(title_str)


#----- body wave

# get arrays of dt(syn-obs), cc 
data_type = "body wave"
data = [ [ float(x[2]), float(x[3]), -1.0*float(x[5]) ] for x in lines if 'surface' not in x[1] ]

weight = np.array([ x[0] for x in data ])
cc = np.array([ x[1] for x in data ])
dt = np.array([ x[2] for x in data ])

idx = np.abs(dt) <= max_dt
dt = dt[idx]
weight = weight[idx]
cc = cc[idx]

# plot CC
nrow = 3
ncol = 1
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of windows')
ax.set_xlim([0.5, 1])
title_str = "%s %.4f(%.2f)" % (data_type, np.sum(cc*weight)/np.sum(weight), np.sum(weight))
ax.set_title(title_str)

# plot dt
nrow = 3
ncol = 2
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(dt, nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of windows')
ax.set_xlim([-max_dt, max_dt])
title_str = "%s %.3f$\pm$%.3f" % (data_type, np.mean(dt), np.std(dt))
ax.set_title(title_str)


#----- All

# get arrays of dt(syn-obs), cc 
data_type = "all"
data = [ [ float(x[2]), float(x[3]), -1.0*float(x[5]) ] for x in lines ]

weight = np.array([ x[0] for x in data ])
cc = np.array([ x[1] for x in data ])
dt = np.array([ x[2] for x in data ])

idx = np.abs(dt) <= max_dt
dt = dt[idx]
weight = weight[idx]
cc = cc[idx]

# plot CC
nrow = 4
ncol = 1
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(cc, nbins, histtype='step')
ax.set_xlabel('cc0')
ax.set_ylabel('No. of windows')
ax.set_xlim([0.5, 1])
title_str = "%s %.4f(%.2f)" % (data_type, np.sum(cc*weight)/np.sum(weight), np.sum(weight))
ax.set_title(title_str)

# plot dt
nrow = 4
ncol = 2
subplot_origin = [ncol-1, nrow-1]*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(dt, nbins, histtype='step')
ax.set_xlabel('dt_cc (s): syn-obs')
ax.set_ylabel('No. of windows')
ax.set_xlim([-max_dt, max_dt])
title_str = "%s %.3f$\pm$%.3f" % (data_type, np.mean(dt), np.std(dt))
ax.set_title(title_str)


#------ save figure
plt.savefig(out_file, format='pdf')