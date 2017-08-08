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
mesh_center_rot = -11.0 # positive value: anti-clockwise rotation viewed from above

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

out_list = "slice_gcircle.list"

R_EARTH_KM = 6371.0
wgs84_sq_one_minus_f = (1.0 - 1.0/298.257223563)**0.5

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
#ndepth = len(depth_list)

#npts = nxi*neta*ndepth
#xyz = np.zeros((3,npts))
#model_value = np.zeros(npts)
#
#dangle = np.deg2rad(dangle)

i = 0
xi = 0.0
xsection_lat = np.zeros(neta)
xsection_lon = np.zeros(neta)
xsection_az = np.zeros(neta)
for eta in eta_list:
  eta = np.deg2rad(eta)
  v1 = np.cos(xi)*vp0 - np.sin(xi)*vr0
  v2 = np.cos(eta)*va0 - np.sin(eta)*vr0
  vx = np.cross(v1,v2)
  vx = vx/sum(vx**2)**0.5
  # local northing, easting unit vectors
  lat = np.arcsin(vx[2])
  lon = np.arctan2(vx[1], vx[0])
  ve = np.array([ -1*np.sin(lon), np.cos(lon), 0.0 ])
  vn = np.array([ -1*np.sin(lat)*np.cos(lon), -1*np.sin(lat)*np.sin(lon), np.cos(lat) ])
  az = np.arctan2(np.dot(v2, ve), np.dot(v2,vn))
  # use geodesic latitude for xsection origin point
  lat_geod = np.arctan2(vx[2], (vx[0]**2+vx[1]**2)**0.5*wgs84_sq_one_minus_f)
  # xsection
  xsection_lat[i] = np.rad2deg(lat_geod)
  xsection_lon[i] = np.rad2deg(lon)
  xsection_az[i] = np.rad2deg(az) + 90.0
  i += 1

#====== write out grid file
with open(out_list, "w") as f:
  for i in range(neta):
    #f.write("%05.2f  %06.2f  %05.1f -25 25 501 5371.00 6371.00 101 EW_lat%05.2f\n" % (xsection_lat[i], xsection_lon[i], xsection_az[i], xsection_lat[i]))
    f.write("%05.2f  %06.2f  %05.1f -20 20 401 5371.00 6371.00 101 EW_lat%05.2f\n" % (xsection_lat[i], xsection_lon[i], xsection_az[i], xsection_lat[i]))