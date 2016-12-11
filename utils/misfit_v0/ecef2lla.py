#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import pyproj

cmt_file = str(sys.argv[1])

# initialize pyproj objects
geod = pyproj.Geod(ellps='WGS84')
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

with open(cmt_file, 'r') as f:
  lines = [ x for x in f.readlines() if not(x.startswith('#')) ]
lines = [x.split(":") for x in lines]
x = float(lines[4][1])
y = float(lines[5][1])
z = float(lines[6][1])

# convert from ECEF(meters) to lla
lon, lat, alt = pyproj.transform(ecef, lla, x, y, z)
dep = -alt / 1000.0

print(lat, " ", lon, " ", dep)