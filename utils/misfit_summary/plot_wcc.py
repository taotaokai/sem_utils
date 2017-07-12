#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

niter = int(sys.argv[1])
misfit_list = sys.argv[2]
out_fig = sys.argv[3]

#====== read misfit list
with open(misfit_list, "r") as f:
  # ignore blank/commeted lines
  lines = [ l for l in (line.strip() for line in f) if l and not(l.startswith('#')) ]
  # split each line into a list 
  lines = [ l.replace('\n','').split() for l in lines ]

nevt = len(lines)
weight = np.zeros((nevt, niter))
wcc = np.zeros((nevt, niter))
wdt = np.zeros((nevt, niter))
evid = [None]*nevt
for ievt in range(nevt):
  l = lines[ievt]
  evid[ievt] = l[0]
  for iter_num in range(niter):
    weight[ievt, iter_num] = float(l[3+iter_num*5])
    wcc[ievt, iter_num] = float(l[4+iter_num*5])
    wdt[ievt, iter_num] = float(l[5+iter_num*5])

#====== plot
iter_nums = range(niter)

fig = plt.figure(figsize=(8.5,11)) # US Letter
ax = fig.add_axes([0.1, 0.05, 0.7, 0.9])

for ievt in range(nevt):
  ax.plot(iter_nums, wcc[ievt,:], 'o-', clip_on=False, lw=0.5, markersize=5, markeredgecolor='none')
  ax.text(niter-1, wcc[ievt,niter-1], "   "+evid[ievt], 
      fontsize=5,
      horizontalalignment='left',
      verticalalignment='center',)
  if wcc[ievt,-1] < np.max(wcc[ievt,:]):
    print("wcc dereases {:s} {:s}: ".format(misfit_list, evid[ievt]), wcc[ievt,:])
    
## total
weight_all = np.sum(weight, axis=0)
wcc_all = np.sum(wcc*weight, axis=0) / weight_all

ax.plot(iter_nums, wcc_all, 'r*', markersize=20, clip_on=False, markerfacecolor='none')

plt.xlabel('No. iteration')
plt.ylabel('wcc')
#str_title = "weighted cc of linearized seismograms: step %.2f wcc_avg %.3f" % (
#    step_length[idx_max], wcc_sum_total[idx_max]/weight_sum_total)
plt.title(misfit_list)

#plt.show()
plt.savefig(out_fig, format='pdf')