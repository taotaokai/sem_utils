#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

fname_list = sys.argv[1]

with open(fname_list, "r") as f:
  # ignore blank/commeted lines
  fnames = [ l for l in (line.strip() for line in f) if l and not(l.startswith('#')) ]

plt.figure(figsize=(11,8.5)) # US Letter

wcc_sum_total = 0.0
weight_sum_total = 0.0

for fname in fnames:

  with open(fname, "r") as f:
    l = f.readline().replace('\n','').split()
    weight_sum = float(l[2])

  with open(fname, "r") as f:
    # ignore blank/commeted lines
    lines = [ l for l in (line.strip() for line in f) if l and not(l.startswith('#')) ]
    # split each line into a list 
    lines = [ l.replace('\n','').split() for l in lines ]

  step_length = np.array([ float(x[0]) for x in lines ])
  wcc_avg = np.array([ float(x[1]) for x in lines ])
  wcc_sum_total += wcc_avg * weight_sum
  weight_sum_total += weight_sum

  plt.plot(step_length, wcc_avg, 'k', lw=0.5)
  idx_max = np.argmax(wcc_avg)
  plt.plot(step_length[idx_max], wcc_avg[idx_max], 'bo', markersize=5)
  plt.text(step_length[idx_max], wcc_avg[idx_max], fname.split('.')[0], fontsize=8)
  if step_length[idx_max] <= 0:
    print("%s step_length[idx_max]=%f" % (fname, step_length[idx_max]))

# total
plt.plot(step_length, wcc_sum_total/weight_sum_total, 'r', lw=2)

idx_max = np.argmax(wcc_sum_total)
plt.plot(step_length[idx_max], wcc_sum_total[idx_max]/weight_sum_total, 'r*', markersize=20)
print("optimal step_length = ", step_length[idx_max])

plt.xlabel('step length')
plt.ylabel('wcc_sum')
str_title = "weighted cc of linearized seismograms: step %.3f wcc_avg %.3f" % (
    step_length[idx_max], wcc_sum_total[idx_max]/weight_sum_total)
plt.title(str_title)


#plt.show()
plt.savefig("wcc_sum.pdf", format='pdf')