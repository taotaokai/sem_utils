#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
from obspy import UTCDateTime

# read command line args
cmt_file = str(sys.argv[1])
srcfrechet_file = str(sys.argv[2])
max_dxs_ratio = float(sys.argv[3])
dcmt_file = str(sys.argv[4])

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
with open(srcfrechet_file, 'r') as f:
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

# project dchi_dmt to be orthogonal with mt
dchi_dmt = dchi_dmt - mt*np.sum(dchi_dmt*mt)/np.sum(mt**2)
cc = np.sum(dchi_dmt*mt)/np.sum(mt**2)**0.5/np.sum(dchi_dmt**2)**0.5
print("cc(dchi_dmt_ortho,mt) = ", cc)

#====== get gradient for xs_ratio and mt_ratio
# xs = R_earth * xs_ratio
# mt = m0 * mt_ratio

dchi_dxs_ratio = R_earth * dchi_dxs

m0 = (0.5*np.sum(mt**2))**0.5
dchi_dmt_ratio = m0 * dchi_dmt

#====== scale CMT gradient
scale_factor = max_dxs_ratio/(np.sum(dchi_dxs_ratio**2))**0.5
dxs_ratio = scale_factor * dchi_dxs_ratio
dmt_ratio = scale_factor * dchi_dmt_ratio

#====== write out new CMTSOLUTION
#xs_perturb = xs + R_earth*dxs_ratio
#
## force mt_perturb to have the same scalar moment as mt 
#mt_perturb = mt + m0*dmt_ratio
#mt_perturb = m0 * mt_perturb/(0.5*np.sum(mt_perturb**2))**0.5

#with open(perturb_cmt_file, 'w') as fp:
#  fp.write('%s\n' % header)
#  fp.write('%-18s %s_perturb\n' % ('event name:', event_id))
#  fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
#  fp.write('%-18s %+15.8E\n' % ('tau(s):',   0.0))
#  fp.write('%-18s %+15.8E\n' % ('x(m):',     xs_perturb[0]))
#  fp.write('%-18s %+15.8E\n' % ('y(m):',     xs_perturb[1]))
#  fp.write('%-18s %+15.8E\n' % ('z(m):',     xs_perturb[2]))
#  fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mt_perturb[0,0]))
#  fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', mt_perturb[1,1]))
#  fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mt_perturb[2,2]))
#  fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mt_perturb[0,1]))
#  fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mt_perturb[0,2]))
#  fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', mt_perturb[1,2]))

# write out dcmt
with open(dcmt_file, 'w') as fp:
  fp.write('%s\n' % header)
  fp.write('%-18s %s_dcmt\n' % ('event name:', event_id))
  fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('tau(s):',   0.0))
  fp.write('%-18s %+15.8E\n' % ('dx(m):',     R_earth*dxs_ratio[0]))
  fp.write('%-18s %+15.8E\n' % ('dy(m):',     R_earth*dxs_ratio[1]))
  fp.write('%-18s %+15.8E\n' % ('dz(m):',     R_earth*dxs_ratio[2]))
  fp.write('%-18s %+15.8E\n' % ('dmxx(N*m):', m0*dmt_ratio[0,0]))
  fp.write('%-18s %+15.8E\n' % ('dmyy(N*m):', m0*dmt_ratio[1,1]))
  fp.write('%-18s %+15.8E\n' % ('dmzz(N*m):', m0*dmt_ratio[2,2]))
  fp.write('%-18s %+15.8E\n' % ('dmxy(N*m):', m0*dmt_ratio[0,1]))
  fp.write('%-18s %+15.8E\n' % ('dmxz(N*m):', m0*dmt_ratio[0,2]))
  fp.write('%-18s %+15.8E\n' % ('dmyz(N*m):', m0*dmt_ratio[1,2]))
