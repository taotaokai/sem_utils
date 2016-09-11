#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert CMTSOLUTION to have ECEF coordinates and tau
"""
import sys
import numpy as np
from obspy import UTCDateTime
import pyproj

cmt_file = str(sys.argv[1])
out_file = str(sys.argv[2])

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

# initialize pyproj objects
geod = pyproj.Geod(ellps='WGS84')
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

tau = float(lines[3][1])/1.628 # mimic triangle with gaussian
lat = float(lines[4][1])
lon = float(lines[5][1])
dep = float(lines[6][1])
# convert from lla to ECEF(meters)
alt = -1000.0 * dep #NOTE ignore local topography
x, y, z = pyproj.transform(lla, ecef, lon, lat, alt)

# centroid time: t0
isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(
    year, month, day, hour, minute, second)
t0 = UTCDateTime(isotime) + time_shift
# modify origin time in header line to have centroid time 
header[1] = str(t0.year)
header[2] = str(t0.month)
header[3] = str(t0.day)
header[4] = str(t0.hour)
header[5] = str(t0.minute)
header[6] = str(t0.second + 1.0e-6*t0.microsecond)

# moment tensor
#1,2,3 -> r,theta,phi
#harvard cmt use dyn*cm for moment tensor, *1e-7 to N*m
m11 = float( lines[7][1])
m22 = float( lines[8][1])
m33 = float( lines[9][1])
m12 = float(lines[10][1])
m13 = float(lines[11][1])
m23 = float(lines[12][1])
mt_rtp = np.array([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]]) * 1.0e-7

# coordinate transformation matrix (r,theta,phi) to (x,y,z) 
r = (x**2 + y**2 + z**2)**0.5
theta = np.arccos(z/r)
phi = np.arctan2(y, x)

sthe = np.sin(theta)
cthe = np.cos(theta)
sphi = np.sin(phi)
cphi = np.cos(phi)

a = np.array(
    [ [ sthe*cphi, cthe*cphi, -1.0*sphi ],
      [ sthe*sphi, cthe*sphi,      cphi ],
      [ cthe     , -1.0*sthe,      0.0  ] ])

mt_xyz = np.dot(np.dot(a, mt_rtp), np.transpose(a))

# write out new CMTSOLUTION_ECEF
with open(out_file, 'w') as fp:
  fp.write('%s\n'            % ' '.join(header))
  fp.write('%-18s %s\n'      % ('event name:',event_id))
  fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('tau(s):',   tau))
  fp.write('%-18s %+15.8E\n' % ('x(m):',     x))
  fp.write('%-18s %+15.8E\n' % ('y(m):',     y))
  fp.write('%-18s %+15.8E\n' % ('z(m):',     z))
  fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mt_xyz[0,0]))
  fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', mt_xyz[1,1]))
  fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mt_xyz[2,2]))
  fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mt_xyz[0,1]))
  fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mt_xyz[0,2]))
  fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', mt_xyz[1,2]))