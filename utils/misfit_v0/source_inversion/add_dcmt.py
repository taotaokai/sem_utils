#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from obspy import UTCDateTime

# read command line args
cmt_file = str(sys.argv[1])
dcmt_file = str(sys.argv[2])
step_length = float(sys.argv[3])
dt0 = float(sys.argv[4])
dtau = float(sys.argv[5])
new_cmt_file = str(sys.argv[6])

##misfit_file = "misfit/misfit.pkl"
#cmt_file = "output_srcfrechet/CMTSOLUTION"
#srcfrechet_file = "output_srcfrechet/src_frechet.000001"
#max_dxs_ratio = 0.001
#
#out_cmt_file = "DATA/CMTSOLUTION.perturb.iter00"
##out_dmodel_file = "DATA/CMTSOLUTION.dmodel.iter00"

# constants
R_earth = 6371000.0

#====== read in original CMTSOLUTION
with open(cmt_file, 'r') as f:
  lines = [ l for l in f.readlines() if not(l.startswith('#')) ]

header = lines[0].split()
year   = header[1]
month  = header[2]
day    = header[3]
hour   = header[4]
minute = header[5]
second = header[6]

lines = [l.split(":") for l in lines]
event_id = lines[1][1].strip()
time_shift = float(lines[2][1])

# centroid time: t0
isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(
    year, month, day, hour, minute, second)
t0 = UTCDateTime(isotime) + time_shift + dt0

# modify origin time in header line to have centroid time 
header[1] = "{:04d}".format(t0.year)
header[2] = "{:02d}".format(t0.month)
header[3] = "{:02d}".format(t0.day)
header[4] = "{:02d}".format(t0.hour)
header[5] = "{:02d}".format(t0.minute)
header[6] = "{:07.4f}".format(t0.second + 1.0e-6*t0.microsecond)
header = ' '.join(header)

# source width
tau = float(lines[3][1]) + dtau

# read xs
x   = float(lines[4][1])
y   = float(lines[5][1])
z   = float(lines[6][1])
# read mt
mxx = float(lines[7][1])
myy = float(lines[8][1])
mzz = float(lines[9][1])
mxy = float(lines[10][1])
mxz = float(lines[11][1])
myz = float(lines[12][1])

#====== read in dcmt
with open(dcmt_file, 'r') as f:
  lines = [ l for l in f.readlines() if not(l.startswith('#')) ]

lines = [l.split(":") for l in lines]
dx   = float(lines[4][1])
dy   = float(lines[5][1])
dz   = float(lines[6][1])
dmxx = float(lines[7][1])
dmyy = float(lines[8][1])
dmzz = float(lines[9][1])
dmxy = float(lines[10][1])
dmxz = float(lines[11][1])
dmyz = float(lines[12][1])

#====== write out new CMTSOLUTION
# force mt_perturb to have the same scalar moment as mt 
#mt_perturb = mt + m0*dmt_ratio
#mt_perturb = m0 * mt_perturb/(0.5*np.sum(mt_perturb**2))**0.5

with open(new_cmt_file, 'w') as fp:
  fp.write('%s\n' % header)
  fp.write('%-18s %s\n' % ('event name:', event_id))
  fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('tau(s):',   tau))
  fp.write('%-18s %+15.8E\n' % ('x(m):',     x + step_length*dx))
  fp.write('%-18s %+15.8E\n' % ('y(m):',     y + step_length*dy))
  fp.write('%-18s %+15.8E\n' % ('z(m):',     z + step_length*dz))
  fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mxx + step_length*dmxx))
  fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', myy + step_length*dmyy))
  fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mzz + step_length*dmzz))
  fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mxy + step_length*dmxy))
  fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mxz + step_length*dmxz))
  fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', myz + step_length*dmyz))
