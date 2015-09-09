#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime

#====== parameters
window_file = str(sys.argv[1])

#====== read window list
with open(window_file, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.replace('\n','').split('|') for x in lines]
window_dicts = [ 
        {'network':     x[0],
         'station':     x[1],
         'location':    x[2],
         'channel':     x[3],
         'azimuth':     float(x[4]),
         'dip':         float(x[5]),
         'starttime':   UTCDateTime(x[6]),
         'endtime':     UTCDateTime(x[7]), 
         'tshift':      float(x[8]),
         'cc_max':      float(x[9]),
         'cc_0':        float(x[10]),
         'ar_0':        float(x[11]),
         } for x in lines ]

stnm   = [ x['station'][2:4]    for x in window_dicts]
tshift = np.array([ x['tshift'] for x in window_dicts])
cc_max = np.array([ x['cc_max'] for x in window_dicts])
cc_0   = np.array([ x['cc_0']   for x in window_dicts])
ar_0   = np.array([ x['ar_0']   for x in window_dicts])

#====== plot tshift/cc_max

idx = cc_max>0.8

x = range(sum(idx))

stnm_idx = [stnm[i] for i in range(len(stnm)) if idx[i] ]

fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

ax0.plot(x, tshift[idx], 'k-*')
ax0.grid()

ax1.plot(x, cc_max[idx], 'k-*', x, cc_0[idx], 'r-o')
ax1.grid()

plt.xticks(x, stnm_idx)
plt.show()

# output results
print 'mean(cc_max) ', np.mean(cc_max[idx])
print 'mean(tshift) ', np.mean(tshift[idx])