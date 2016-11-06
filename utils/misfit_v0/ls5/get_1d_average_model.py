#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np 

import matplotlib
import matplotlib.pyplot as plt

# input
data_file = str(sys.argv[1])

# parameters
min_dr = 2.0e-5
nmodel = 4

# read data file
with open(data_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]

r   = np.array([float(l[0]) for l in lines])
r_center = np.array([float(l[1]) for l in lines])
nr = r.size

model = np.zeros((nr, nmodel))
for i in range(nmodel):
  model[:,i] = np.array([float(l[2+i]) for l in lines])

# sort radius
ix = np.argsort(r)
r = r[ix]
r_center = r_center[ix]
model = model[ix,:]

# group points with almost the same radius
# and take the mean value
dr = np.diff(r)
indices = np.where(dr > min_dr)
indices = indices[0]

n = indices.size
r_1d = np.zeros(n+2)
model_1d = np.zeros((n+2, nmodel))

idx1 = 0
idx_pad = 0
for i in range(n):
  idx2 = indices[i] + 1

  r1 = r[idx1:idx2]
  r1_center = r_center[idx1:idx2]
  model1 = model[idx1:idx2,:]

  # separate 410-km
  r1_mean = np.mean(r1)
  if np.abs(r1_mean - 0.93564) < min_dr:
    # below
    ix2 = r1_center < 0.93564
    r_1d[i+idx_pad] = np.mean(r1[ix2])
    model_1d[i+idx_pad] = np.mean(model1[ix2,:], axis=0)

    idx_pad += 1

    # above
    ix2 = r1_center > 0.93564
    r_1d[i+idx_pad] = np.mean(r1[ix2])
    model_1d[i+idx_pad] = np.mean(model1[ix2,:], axis=0)

  # separate 660-km
  elif np.abs(r1_mean - 0.89797) < min_dr:
    # below
    ix2 = r1_center < 0.89797
    r_1d[i+idx_pad] = np.mean(r1[ix2])
    model_1d[i+idx_pad] = np.mean(model1[ix2,:], axis=0)

    idx_pad += 1

    # above
    ix2 = r1_center > 0.89797
    r_1d[i+idx_pad] = np.mean(r1[ix2])
    model_1d[i+idx_pad] = np.mean(model1[ix2,:], axis=0)

  else:
    r_1d[i+idx_pad] = np.mean(r1)
    model_1d[i+idx_pad,:] = np.mean(model1, axis=0)

  idx1 = idx2

# write out average 1d profile
out_file = "%s_1d.txt" % (data_file.split(".")[0])
with open(out_file, 'w') as f:
  for i in range(n+2):
    f.write("%18.8e  " % (r_1d[i]))
    for imodel in range(nmodel):
      f.write("%18.8e  " % (model_1d[i,imodel]) )
    f.write("\n")

# plot
#plt.figure()
#plt.plot(model_1d[:,0], r_1d, 'b.-')
#plt.show()