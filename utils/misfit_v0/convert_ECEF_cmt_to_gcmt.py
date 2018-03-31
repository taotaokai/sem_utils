#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert CMTSOLUTION from ECEF format to gcmt format
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
GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

tau = float(lines[3][1])
x   = float(lines[4][1])
y   = float(lines[5][1])
z   = float(lines[6][1])
# convert from ECEF(meters) to lla
lon, lat, alt = pyproj.transform(ecef, lla, x, y, z)
dep = -alt / 1000.0

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

# moment tensor
# ECEF=false: 1,2,3 -> r,theta,phi
# ECEF=true:  1,2,3 -> x,y,z
m11 = float( lines[7][1])
m22 = float( lines[8][1])
m33 = float( lines[9][1])
m12 = float(lines[10][1])
m13 = float(lines[11][1])
m23 = float(lines[12][1])
mt = np.array([
  [m11, m12, m13], 
  [m12, m22, m23], 
  [m13, m23, m33]])
# transform from spherical to cartesian coordinate
r = (x**2 + y**2 + z**2)**0.5
theta = np.arccos(z/r)
phi = np.arctan2(y, x)
# rotation matrix
sthe = np.sin(theta)
cthe = np.cos(theta)
sphi = np.sin(phi)
cphi = np.cos(phi)
# basis transform matrix: e_x,y,z = a * e_r,t,p
mt_rtp = np.zeros((3,3))
a = np.array(
    [ [ sthe*cphi, cthe*cphi, -1.0*sphi ],
      [ sthe*sphi, cthe*sphi,      cphi ],
      [ cthe     , -1.0*sthe,      0.0  ] ])
mt_rtp = np.dot(np.dot(np.transpose(a), mt), a)

with open(out_file, 'w') as fp:
  fp.write('%s\n'            % ' '.join(header))
  fp.write('%-18s %s\n'      % ('event name:', event_id))
  fp.write('%-18s %+15.8E\n' % ('time shift:',    0.0))
  fp.write('%-18s %+15.8E\n' % ('half duration:', tau*1.628))
  fp.write('%-18s %+15.8E\n' % ('latitude:',    lat))
  fp.write('%-18s %+15.8E\n' % ('longitude:',   lon))
  fp.write('%-18s %+15.8E\n' % ('depth:',       dep))
  fp.write('%-18s %+15.8E\n' % ('Mrr(dyn*cm):', mt_rtp[0,0]*1.0e7))
  fp.write('%-18s %+15.8E\n' % ('Mtt(dyn*cm):', mt_rtp[1,1]*1.0e7))
  fp.write('%-18s %+15.8E\n' % ('Mpp(dyn*cm):', mt_rtp[2,2]*1.0e7))
  fp.write('%-18s %+15.8E\n' % ('Mrt(dyn*cm):', mt_rtp[0,1]*1.0e7))
  fp.write('%-18s %+15.8E\n' % ('Mrp(dyn*cm):', mt_rtp[0,2]*1.0e7))
  fp.write('%-18s %+15.8E\n' % ('Mtp(dyn*cm):', mt_rtp[1,2]*1.0e7))
