#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert CMTSOLUTION in ECEF coordinate to Geodetic coord.
"""
import sys
import numpy as np
# from scipy.interpolate import interpn

# from netCDF4 import Dataset
from obspy import UTCDateTime

from pyproj import Transformer

#====== user input
cmt_file = str(sys.argv[1])
out_file = str(sys.argv[2])
# topoGRDfile = str(sys.argv[3]) # 'ETOPO1_Ice_g_smooth_b15km_I4m.grd'

epsg_ecef = "EPSG:4978"
epsg_wgs84 = "EPSG:4326"
transformer = Transformer.from_crs(epsg_ecef, epsg_wgs84)

#--- read in CMTSOLUTION in ECEF coord.
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

hdur_triangle = float(lines[3][1])*1.628 # mimic triangle with gaussian
x = float(lines[4][1]) # meter
y = float(lines[5][1])
z = float(lines[6][1])

# convert from ECEF to llh
lat, lon, height = transformer.transform(x, y, z)
# topo_local = interpn((grd_lon,grd_lat), grd_topo, [lon,lat])
# print("lon,lat,topo = ",lon,lat,topo_local)
# dep_km = (topo_local - height)/1000.0 # elliptic height, ellipsoid approx. MSL
dep_km = (0 - height)/1000.0 # elliptic height, ellipsoid approx. MSL

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
#1,2,3 -> x,y,z
#harvard cmt use dyn*cm for moment tensor, 1 N*m = 1.0e7 dyn*cm
m11 = float( lines[7][1])
m22 = float( lines[8][1])
m33 = float( lines[9][1])
m12 = float(lines[10][1])
m13 = float(lines[11][1])
m23 = float(lines[12][1])
mt_xyz = np.array([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]]) * 1.0e7

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

mt_rtp = np.dot(np.dot(np.transpose(a), mt_xyz), a)

# write out new CMTSOLUTION_ECEF
with open(out_file, 'w') as fp:
  fp.write('%s\n'            % ' '.join(header))
  fp.write('%-18s %s\n'      % ('event_name:',       event_id))
  fp.write('%-18s %+15.8E\n' % ('time_shift(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('hdur_triangle(s):', hdur_triangle))
  fp.write('%-18s %+15.8E\n' % ('latitude(deg):',    lat))
  fp.write('%-18s %+15.8E\n' % ('longitude(deg):',   lon))
  fp.write('%-18s %+15.8E\n' % ('depth(km):',        dep_km))
  fp.write('%-18s %+15.8E\n' % ('Mrr(dyn*cm):', mt_rtp[0,0]))
  fp.write('%-18s %+15.8E\n' % ('Mtt(dyn*cm):', mt_rtp[1,1]))
  fp.write('%-18s %+15.8E\n' % ('Mpp(dyn*cm):', mt_rtp[2,2]))
  fp.write('%-18s %+15.8E\n' % ('Mrt(dyn*cm):', mt_rtp[0,1]))
  fp.write('%-18s %+15.8E\n' % ('Mrp(dyn*cm):', mt_rtp[0,2]))
  fp.write('%-18s %+15.8E\n' % ('Mtp(dyn*cm):', mt_rtp[1,2]))
