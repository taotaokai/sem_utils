#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Plot .nc files generated by xsem_slice_gcircle
Need both inverted and reference model
"""
import sys
import numpy as np

from netCDF4 import Dataset
import pyproj

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#------ user inputs
nc_file_list = sys.argv[1]
title = sys.argv[2]
out_fig = sys.argv[3]

marker_interval = 5.0 #degree
#marker_theta = np.arange(-1,5)*np.deg2rad(marker_interval) # Ryukyu
marker_theta = np.arange(0,5)*np.deg2rad(marker_interval) # Japan

# earth parameters 
R_earth_meter = 6371000.0
R_earth_km = 6371.0

#-- plot controls
# map plot
map_lat_center = 38.5
map_lon_center = 118.0
map_lat_min=  5
map_lat_max= 55
map_lon_min= 90
map_lon_max= 160
map_parallels = np.arange(0.,81,10.)
map_meridians = np.arange(0.,351,10.)

# model plot region
width = 0.8
height = 0.6

#------ utility functions
def rotation_matrix(v_axis, theta):
  """ rotation matrix: rotate through a given axis by an angle of theta
      (right-hand rule)
  """
  sint = np.sin(theta)
  cost = np.cos(theta)
  one_minus_cost = 1.0 - cost

  # normalize v_axis to unit vector
  v = v_axis / sum(v_axis**2)**0.5

  # rotation matrix
  R = np.zeros((3,3))
  R[0,0] = cost + one_minus_cost * v_axis[0]**2
  R[1,1] = cost + one_minus_cost * v_axis[1]**2
  R[2,2] = cost + one_minus_cost * v_axis[2]**2

  R[0,1] = one_minus_cost*v_axis[0]*v_axis[1] - sint * v_axis[2]
  R[1,0] = one_minus_cost*v_axis[0]*v_axis[1] + sint * v_axis[2]

  R[0,2] = one_minus_cost*v_axis[0]*v_axis[2] + sint * v_axis[1]
  R[2,0] = one_minus_cost*v_axis[0]*v_axis[2] - sint * v_axis[1]

  R[1,2] = one_minus_cost*v_axis[1]*v_axis[2] - sint * v_axis[0]
  R[2,1] = one_minus_cost*v_axis[1]*v_axis[2] + sint * v_axis[0]

  return R

#------ plot map and xsection surface trace and marker
fig = plt.figure(figsize=(8.5,11))

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
m = Basemap(ax=ax, projection='tmerc', 
    resolution='l', area_thresh=30000,
    llcrnrlat=map_lat_min, llcrnrlon=map_lon_min, 
    urcrnrlat=map_lat_max, urcrnrlon=map_lon_max,
    lat_0=map_lat_center, lon_0=map_lon_center,
    )

m.drawcoastlines(linewidth=0.2)
m.drawcountries(linewidth=0.2)
m.drawparallels(map_parallels, linewidth=0.1, labels=[1,0,0,0], fontsize=24)
m.drawmeridians(map_meridians, linewidth=0.1, labels=[0,0,0,1], fontsize=24)
#m.shadedrelief()

# initialize pyproj objects
geod = pyproj.Geod(ellps='WGS84')
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

#--- plot Holocene volcanoes
volcano_file= 'volcanoes.list'
with open(volcano_file, 'r') as f:
  lines = [ l.split('|') for l in f.readlines() if not l.startswith('#') ]
  volcano_lats = np.array([float(x[4]) for x in lines])
  volcano_lons = np.array([float(x[5]) for x in lines])

x, y = m(volcano_lons, volcano_lats)
ax.plot(x, y, '^', color='none', markersize=14, markeredgecolor='red')

#-- plot geological blocks 
block_line_file = 'zhangpz_block.txt'
block_lines = []
with open(block_line_file, 'r') as f:
  lon = []
  lat = []
  for l in f.readlines():
    if not l.startswith('>'):
      x = l.split()
      lon.append(float(x[0]))
      lat.append(float(x[1]))
    else:
      block_lines.append([lon, lat])
      lon = []
      lat = []
for l in block_lines:
  x, y = m(l[0], l[1])
  ax.plot(x, y, lw=0.5, color='gray')

#-- plot plate_boundary
pb_line_file = 'zhangpz_pb.txt'
pb_lines = []
with open(pb_line_file, 'r') as f:
  lon = []
  lat = []
  for l in f.readlines():
    if not l.startswith('>'):
      x = l.split()
      lon.append(float(x[0]))
      lat.append(float(x[1]))
    else:
      pb_lines.append([lon, lat])
      lon = []
      lat = []
for l in pb_lines:
  x, y = m(l[0], l[1])
  ax.plot(x, y, lw=1.0, color='red')

#------ plot xsection surfae tracks 

#--- read nc list
with open(nc_file_list, 'r') as f:
  nc_list = [ l.split()[0] for l in f.readlines() if not l.startswith('#') ]

for nc_file in nc_list:
  print(nc_file)

  # read nc files 
  fh = Dataset(nc_file, mode='r')
  
  radius = fh.variables['radius'][:]
  theta = np.deg2rad(fh.variables['theta'][:])
  
  xsection_lat0 = fh.variables['latitude_orig'][:]
  xsection_lon0 = fh.variables['longitude_orig'][:]
  xsection_lat_rot = fh.variables['latitude_rotation'][:]
  xsection_lon_rot = fh.variables['longitude_rotation'][:]
  xsection_ang_rot = np.deg2rad(fh.variables['rotation_angle'][:])
  xsection_lats = fh.variables['latitude'][:]
  xsection_lons = fh.variables['longitude'][:]

  # xsection geometry parameters
  # unit direction vector v_orig at origin point of the great circle
  x, y, z = pyproj.transform(lla, ecef, xsection_lon0, xsection_lat0, 0.0)
  v_orig = np.array([x, y, z])
  v_orig = v_orig / np.sqrt(np.sum(v_orig**2))
  
  # rotation axis
  x, y, z = pyproj.transform(lla, ecef, xsection_lon_rot, xsection_lat_rot, 0.0)
  v_axis = np.array([x, y, z])
  v_axis = v_axis / np.sqrt(np.sum(v_axis**2))
  
  cone_angle = np.arccos(np.dot(v_orig, v_axis))

  # plot xsection surface track
  x, y = m(xsection_lons, xsection_lats)
  ax.plot(x, y, 'k-', lw=0.5)
  xlim = ax.get_xlim()
  length = xlim[1]-xlim[0]
  ax.arrow(x[-1], y[-1], x[-1]-x[-2], y[-1]-y[-2], head_width=0.02*length, head_length=0.03*length, fc='k', ec='k')

  # plot xsection surface marker 
  nmarker = len(marker_theta)
  marker_lons = np.zeros(nmarker)
  marker_lats = np.zeros(nmarker)
  for i in range(nmarker):
    rotmat = rotation_matrix(v_axis, marker_theta[i]/np.sin(cone_angle))
    vr = np.dot(rotmat, v_orig)*R_earth_meter
    marker_lons[i], marker_lats[i], alt = pyproj.transform(ecef, lla, vr[0], vr[1], vr[2])
  x, y = m(marker_lons, marker_lats)
  ax.plot(x, y, 'ro', markersize=8, )

# title
ax.set_title(title)

#------ save figure
#plt.show()
plt.savefig(out_fig, format='pdf')