#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
#
import numpy as np
#
import pyproj
#
#import matplotlib
##matplotlib.use("pdf")
#import matplotlib.pyplot as plt

#====== paramters
mesh_center_lat = 35.2
mesh_center_lon = 116.8
mesh_center_rot = -11.0 # anti-clockwise rotation viewed from above

# when zero rotation angle: xi -> easting, eta -> northing
xi0 = -30
xi1 = 30
eta0 = -30
eta1 = 30
dangle = 5.0

depth0 = 100
depth1 = 900
ddepth = 200.0

peak_value = 0.01

out_list = "xyz.list"

R_EARTH_KM = 6371.0

#====== get ECEF coordinate for the center point
# initialize pyproj objects
geod = pyproj.Geod(ellps='WGS84')
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

# local radial unit vector
x, y, z = pyproj.transform(lla, ecef, mesh_center_lon, mesh_center_lat, 0.0)
vr0 = np.array([x,y,z])
vr0 /= sum(vr0**2)**0.5

# local northing, easting unit vectors
lat0 = np.arcsin(vr0[2])
lon0 = np.arctan2(vr0[1], vr0[0])
ve0 = np.array([ -1*np.sin(lon0), np.cos(lon0), 0.0 ])
vn0 = np.array([ -1*np.sin(lat0)*np.cos(lon0), -1*np.sin(lat0)*np.sin(lon0), np.cos(lat0) ])

# local basis along and perpendicular to the rotation azimuth
az0 = -1*np.deg2rad(mesh_center_rot)
va0 =    vn0*np.cos(az0) + ve0*np.sin(az0)
vp0 = -1*vn0*np.sin(az0) + ve0*np.cos(az0)

#====== get grid of xi, eta: xi along va0, eta along vp0
xi_list = np.arange(xi0,xi1+dangle/2,dangle)
eta_list = np.arange(eta0,eta1+dangle/2,dangle)
depth_list = np.arange(depth0,depth1+ddepth/2,ddepth)

nxi = len(xi_list)
neta = len(eta_list)
ndepth = len(depth_list)

npts = nxi*neta*ndepth
xyz = np.zeros((3,npts))
model_value = np.zeros(npts)

dangle = np.deg2rad(dangle)

i = 0
for xi in xi_list:
  xi = np.deg2rad(xi)
  for eta in eta_list:
    eta = np.deg2rad(eta)
    v1 = np.cos(xi)*vp0 - np.sin(xi)*vr0
    v2 = np.cos(eta)*va0 - np.sin(eta)*vr0
    vx = np.cross(v1,v2)
    vx = vx/sum(vx**2)**0.5
    for depth in depth_list:
      sign = np.cos(np.pi*xi/dangle)*np.cos(np.pi*eta/dangle)*np.cos(np.pi*(depth-depth0)/ddepth)
      r = 1.0 - depth/R_EARTH_KM
      xyz[:,i] = r*vx
      model_value[i] = sign*peak_value
      i += 1

#====== write out grid file
with open(out_list, "w") as f:
  for i in range(npts):
    f.write("%+14.7E  %+14.7E  %+14.7E  40  %+8.5f\n" % (xyz[0,i], xyz[1,i], xyz[2,i], model_value[i]))