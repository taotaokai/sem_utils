#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pyproj


GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

# slicev2

lon = 104.484203
lat =  32.497388
alt = 0.0
x1, y1, z1 = pyproj.transform(lla, ecef, lon, lat, alt)

lon = 130.116009
lat =  22.687918
alt = 0.0
x2, y2, z2 = pyproj.transform(lla, ecef, lon, lat, alt)

v1 = np.array([x1, y1, z1])
v2 = np.array([x2, y2, z2])
vn = np.cross(v1, v2)

vn = vn/sum(vn**2)**0.5

print(vn)

## west edge
#
#lon = 127.239369
#lat =  54.882688
#alt = 0.0
#x1, y1, z1 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#lon = 112.676065
#lat =  30.476822
#alt = 0.0
#x2, y2, z2 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#v1 = np.array([x1, y1, z1])
#v2 = np.array([x2, y2, z2])
#vn = np.cross(v1, v2)
#
#vn = vn/sum(vn**2)**0.5
#
#print(vn)
#
## south edge
#
#lon = 111.576935
#lat =  27.163477
#alt = 0.0
#x1, y1, z1 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#lon = 134.635663
#lat =  16.973191
#alt = 0.0
#x2, y2, z2 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#v1 = np.array([x1, y1, z1])
#v2 = np.array([x2, y2, z2])
#vn = np.cross(v1, v2)
#
#vn = vn/sum(vn**2)**0.5
#
#print(vn)
#
## north edge
#
#lon = 153.107141
#lat =  40.186350
#alt = 0.0
#x1, y1, z1 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#lon = 127.239369
#lat =  54.882688
#alt = 0.0
#x2, y2, z2 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#v1 = np.array([x1, y1, z1])
#v2 = np.array([x2, y2, z2])
#vn = np.cross(v1, v2)
#
#vn = vn/sum(vn**2)**0.5
#
#print(vn)
#
## east edge
#
#lon = 132.411306
#lat =  18.179365
#alt = 0.0
#x1, y1, z1 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#lon = 152.215893
#lat =  41.126297
#alt = 0.0
#x2, y2, z2 = pyproj.transform(lla, ecef, lon, lat, alt)
#
#v1 = np.array([x1, y1, z1])
#v2 = np.array([x2, y2, z2])
#vn = np.cross(v1, v2)
#
#vn = vn/sum(vn**2)**0.5
#
#print(vn)