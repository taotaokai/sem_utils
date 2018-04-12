#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pyproj

lons = np.arange(90, 150.1, 0.25)
lats = np.arange(10, 60.1, 0.25)
deps = np.arange(0, 901, 10) # km

R_EARTH = 6371000.0

GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

lon_grid, lat_grid, dep_grid = np.meshgrid(lons, lats, deps, indexing='ij')

xx, yy, zz = pyproj.transform(lla, ecef, lon_grid, lat_grid, -1000.0*dep_grid)
xx = xx/R_EARTH
yy = yy/R_EARTH
zz = zz/R_EARTH

#rr = 1.0 - dep_grid/R_earth_km
#xx = np.cos(np.deg2rad(lat_grid)) * np.cos(np.deg2rad(lon_grid)) * rr
#yy = np.cos(np.deg2rad(lat_grid)) * np.sin(np.deg2rad(lon_grid)) * rr
#zz = np.sin(np.deg2rad(lat_grid)) * rr

n = xx.size
aa = np.zeros((n, 3))
aa[:,0] = xx.flatten()
aa[:,1] = yy.flatten()
aa[:,2] = zz.flatten()

np.savetxt('xyz.lst', aa, fmt="%f %f %f")

bb = np.zeros((n, 3))
bb[:,0] = lon_grid.flatten()
bb[:,1] = lat_grid.flatten()
bb[:,2] = dep_grid.flatten()

np.savetxt('lld.lst', bb, fmt="%e %e %e")