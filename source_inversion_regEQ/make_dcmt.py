#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import numpy as np
from obspy import UTCDateTime

parser = argparse.ArgumentParser()

parser.add_argument("cmt_file" )  # "CMTSOLUTION"
parser.add_argument("srcfrechet_file") # srcfrechet.00001
parser.add_argument("max_dxs_ratio", type=float) # 0.001
parser.add_argument("max_dmt_ratio", type=float) # 0.01
parser.add_argument("out_cmt_file_dxs") # dxs.cmt
parser.add_argument("out_cmt_file_dmt") # dmt.cmt

args = parser.parse_args()

#====== read in original CMTSOLUTION
with open(args.cmt_file, 'r') as f:
  lines = [ x for x in f.readlines() if not(x.startswith('#')) ]

header = lines[0].split()
year   = header[1]
month  = header[2]
day    = header[3]
hour   = header[4]
minute = header[5]
second = header[6]

lines = [x.split(":") for x in lines]
event_id = lines[1][1].strip()
time_shift = float(lines[2][1])

# centroid time: t0
isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(
    year, month, day, hour, minute, second)
t0 = UTCDateTime(isotime) + time_shift
# modify origin time in header line to have centroid time
header[1] = "{:04d}".format(t0.year)
header[2] = "{:02d}".format(t0.month)
header[3] = "{:02d}".format(t0.day)
header[4] = "{:02d}".format(t0.hour)
header[5] = "{:02d}".format(t0.minute)
header[6] = "{:07.4f}".format(t0.second + 1.0e-6*t0.microsecond)

header = ' '.join(header)

#====== read in source gradient
with open(args.srcfrechet_file, 'r') as f:
  lines = [ x for x in f.readlines() if not(x.startswith('#')) ]

lines = [x.split() for x in lines]

t0  = float(lines[0][0]);  dchi_dt0  = float(lines[0][1])
tau = float(lines[1][0]);  dchi_dtau = float(lines[1][1])

xs = np.zeros(3); dchi_dxs = np.zeros(3)
xs[0] = float(lines[2][0]);  dchi_dxs[0] = float(lines[2][1])
xs[1] = float(lines[3][0]);  dchi_dxs[1] = float(lines[3][1])
xs[2] = float(lines[4][0]);  dchi_dxs[2] = float(lines[4][1])

mt = np.zeros((3,3)); dchi_dmt = np.zeros((3,3));
mt[0,0] = float(lines[5][0]);  dchi_dmt[0,0] = float(lines[5][1])
mt[1,1] = float(lines[6][0]);  dchi_dmt[1,1] = float(lines[6][1])
mt[2,2] = float(lines[7][0]);  dchi_dmt[2,2] = float(lines[7][1])

mt[0,1] = float(lines[8][0]);  dchi_dmt[0,1] = float(lines[8][1])
mt[1,0] = mt[0,1]; dchi_dmt[1,0] = dchi_dmt[0,1]

mt[0,2] = float(lines[9][0]);  dchi_dmt[0,2] = float(lines[9][1])
mt[2,0] = mt[0,2]; dchi_dmt[2,0] = dchi_dmt[0,2]

mt[1,2] = float(lines[10][0]); dchi_dmt[1,2] = float(lines[10][1])
mt[2,1] = mt[1,2]; dchi_dmt[2,1] = dchi_dmt[1,2]

# correlation between mt and dchi_dmt
cc = np.sum(dchi_dmt*mt)/np.sum(mt**2)**0.5/np.sum(dchi_dmt**2)**0.5
print("cc(dchi_dmt,mt) = ", cc)

## [2017-05-28] I decide not to project dchi_dmt orthogonal to mt
## project dchi_dmt to be orthogonal with mt
#dchi_dmt = dchi_dmt - mt*np.sum(dchi_dmt*mt)/np.sum(mt**2)
#cc = np.sum(dchi_dmt*mt)/np.sum(mt**2)**0.5/np.sum(dchi_dmt**2)**0.5
#print("cc(dchi_dmt_ortho,mt) = ", cc)

#====== scale gradients

r0 = sum(xs**2)**0.5
dxs_scaled = max_dxs_ratio * r0 / sum(dchi_dxs**2)**0.5 * dchi_dxs

m0 = (0.5*sum(mt**2))**0.5
dmt_scaled = max_dmt_ratio * m0 / (0.5*sum(dchi_dmt**2))**0.5 * dchi_dmt

#====== write out new CMTSOLUTION

# write out dcmt
with open(args.out_cmtfile_dxs, 'w') as fp:
  fp.write('%s\n' % header)
  fp.write('%-18s %s_dcmt\n' % ('event_name:', event_id))
  fp.write('%-18s %+15.8E\n' % ('dt0(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('dtau(s):',   0.0))
  fp.write('%-18s %+15.8E\n' % ('dx(m):',     dxs_scaled[0]))
  fp.write('%-18s %+15.8E\n' % ('dy(m):',     dxs_scaled[1]))
  fp.write('%-18s %+15.8E\n' % ('dz(m):',     dxs_scaled[2]))
  fp.write('%-18s %+15.8E\n' % ('dmxx(N*m):', 0.0))
  fp.write('%-18s %+15.8E\n' % ('dmyy(N*m):', 0.0))
  fp.write('%-18s %+15.8E\n' % ('dmzz(N*m):', 0.0))
  fp.write('%-18s %+15.8E\n' % ('dmxy(N*m):', 0.0))
  fp.write('%-18s %+15.8E\n' % ('dmxz(N*m):', 0.0))
  fp.write('%-18s %+15.8E\n' % ('dmyz(N*m):', 0.0))

with open(args.out_cmtfile_dmt, 'w') as fp:
  fp.write('%s\n' % header)
  fp.write('%-18s %s_dcmt\n' % ('event_name:', event_id))
  fp.write('%-18s %+15.8E\n' % ('dt0(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('dtau(s):',   0.0))
  fp.write('%-18s %+15.8E\n' % ('dx(m):',     0.0))
  fp.write('%-18s %+15.8E\n' % ('dy(m):',     0.0))
  fp.write('%-18s %+15.8E\n' % ('dz(m):',     0.0))
  fp.write('%-18s %+15.8E\n' % ('dmxx(N*m):', dmt_scaled[0,0]))
  fp.write('%-18s %+15.8E\n' % ('dmyy(N*m):', dmt_scaled[1,1]))
  fp.write('%-18s %+15.8E\n' % ('dmzz(N*m):', dmt_scaled[2,2]))
  fp.write('%-18s %+15.8E\n' % ('dmxy(N*m):', dmt_scaled[0,1]))
  fp.write('%-18s %+15.8E\n' % ('dmxz(N*m):', dmt_scaled[0,2]))
  fp.write('%-18s %+15.8E\n' % ('dmyz(N*m):', dmt_scaled[1,2]))
