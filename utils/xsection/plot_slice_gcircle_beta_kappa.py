#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" contourfill plot dlnvs with 1D Vp/Vs profiles overlaid
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
nc_file = sys.argv[1]
title = sys.argv[2]
out_fig = sys.argv[3]

model_names = ['vp0', 'vs0', 'beta', 'alpha']

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

# model colormap
#cmap = plt.cm.get_cmap("jet_r")
#cmap.set_under("white")
#cmap.set_over("white")

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

#------ read nc files 
fh = Dataset(nc_file, mode='r')

radius = fh.variables['radius'][:]
theta = np.deg2rad(fh.variables['theta'][:])

xsection_lat0 = fh.variables['latitude_orig'][:]
xsection_lon0 = fh.variables['longitude_orig'][:]
xsection_az0 = fh.variables['azimuth_orig'][:]
xsection_lats = fh.variables['latitude'][:]
xsection_lons = fh.variables['longitude'][:]

model = {}
for tag in model_names:
  model[tag] = fh.variables[tag][:]

vp = (1.0 + model['alpha']) * model['vp0']
vs = (1.0 + model['beta']) * model['vs0']
kappa = vp/vs

#------ plot map and xsection surface trace and marker
fig = plt.figure(figsize=(8.5,11))

ax = fig.add_axes([0.2, 0.75, 0.6, 0.2])
m = Basemap(ax=ax, projection='tmerc', resolution='l',
    llcrnrlat=map_lat_min, llcrnrlon=map_lon_min, 
    urcrnrlat=map_lat_max, urcrnrlon=map_lon_max,
    lat_0=map_lat_center, lon_0=map_lon_center)
m.drawcoastlines(linewidth=0.2)
m.drawcountries(linewidth=0.2)
m.drawparallels(map_parallels, linewidth=0.1, labels=[1,0,0,0], fontsize=8)
m.drawmeridians(map_meridians, linewidth=0.1, labels=[0,0,0,1], fontsize=8)

# initialize pyproj objects
geod = pyproj.Geod(ellps='WGS84')
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

#--- xsection geometry parameters
# unit direction vector v0 at origin point of the great circle
x, y, z = pyproj.transform(lla, ecef, xsection_lon0, xsection_lat0, 0.0)
v0 = np.array([x, y, z])
v0 = v0 / np.sqrt(np.sum(v0**2))

# unit direction vector v1 along the shooting azimuth of the great circle
vnorth = np.array( [ - np.sin(np.deg2rad(xsection_lat0)) * np.cos(np.deg2rad(xsection_lon0)),
                     - np.sin(np.deg2rad(xsection_lat0)) * np.sin(np.deg2rad(xsection_lon0)),
                       np.cos(np.deg2rad(xsection_lat0)) ])
veast = np.array([ - np.sin(np.deg2rad(xsection_lon0)), np.cos(np.deg2rad(xsection_lon0)), 0.0 ])

v1 = np.cos(np.deg2rad(xsection_az0)) * vnorth + np.sin(np.deg2rad(xsection_az0)) * veast
v1 = v1 / np.sqrt(np.sum(v1**2))

# rotation axis = v0 cross-product v1
v_axis = np.cross(v0, v1)
v_axis = v_axis / np.sqrt(np.sum(v_axis**2))

#--- plot fault lines
fault_line_file = 'fault_lines.txt'
fault_lines = []
with open(fault_line_file, 'r') as f:
  lon = []
  lat = []
  for l in f.readlines():
    if not l.startswith('>'):
      x = l.split()
      lon.append(float(x[0]))
      lat.append(float(x[1]))
    else:
      fault_lines.append([lon, lat])
      lon = []
      lat = []
for l in fault_lines:
  x, y = m(l[0], l[1])
  ax.plot(x, y, 'k-', lw=0.1)

#--- plot seismicity
catalog_file = 'isc_d50km.txt'
with open(catalog_file, 'r') as f:
  lines = [ l.split('|') for l in f.readlines() if not l.startswith('#') ] 
  eq_lats = np.array([float(x[2]) for x in lines])
  eq_lons = np.array([float(x[3]) for x in lines])
  eq_deps = 1000.0 * np.array([float(x[4]) for x in lines])
# limit earthquake depth
idx = eq_deps < 410
eq_lats = eq_lats[idx]
eq_lons = eq_lons[idx]
eq_deps = eq_deps[idx]

eq_x, eq_y, eq_z = pyproj.transform(lla, ecef, eq_lons, eq_lats, -1.0*eq_deps, radians=False)

# get earthquakes located within a certain distance to the xsection
eq_dist = v_axis[0]*eq_x + v_axis[1]*eq_y + v_axis[2]*eq_z
eq_indx = np.abs(eq_dist) <= 25000.0
# project to xsection
nsrc = np.sum(eq_indx)
eq_xyz = np.zeros((nsrc,3))
eq_xyz[:,0] = eq_x[eq_indx]
eq_xyz[:,1] = eq_y[eq_indx]
eq_xyz[:,2] = eq_z[eq_indx]
eq_xyz = eq_xyz - np.sum(eq_xyz*v_axis, axis=1, keepdims=True)*v_axis
eq_theta = np.arctan2(np.sum(eq_xyz*v1, axis=1), np.sum(eq_xyz*v0, axis=1))
eq_radius = np.sum(eq_xyz**2, axis=1)**0.5/1000.0

#--- plot xsection surface track
x, y = m(xsection_lons, xsection_lats)
ax.plot(x, y, 'k-', lw=0.5)
xlim = ax.get_xlim()
length = xlim[1]-xlim[0]
ax.arrow(x[-1], y[-1], x[-1]-x[-2], y[-1]-y[-2], head_width=0.03*length, head_length=0.03*length, fc='k', ec='k')

#--- plot xsection surface marker 
nmarker = 5
marker_lons = np.zeros(nmarker)
marker_lats = np.zeros(nmarker)
marker_theta = np.linspace(theta[0], theta[-1], nmarker)
for i in range(nmarker):
  rotmat = rotation_matrix(v_axis, marker_theta[i])
  vr = np.dot(rotmat, v0)*R_earth_meter
  marker_lons[i], marker_lats[i], alt = pyproj.transform(ecef, lla, vr[0], vr[1], vr[2])
x, y = m(marker_lons, marker_lats)
ax.plot(x, y, 'ro', markersize=4, )

# title
ax.set_title(title)

#------ plot models
# shift theta to center plot 
theta_mid = (theta[-1]+theta[0])/2.0
theta = theta - theta_mid

rr, tt = np.meshgrid(radius, theta, indexing='ij')
xx = np.sin(tt) * rr
yy = np.cos(tt) * rr

# seismicity
eq_theta = eq_theta - theta_mid  
idx = (eq_theta >= np.min(theta)) & (eq_theta <= np.max(theta))
eq_x = eq_radius[idx] * np.sin(eq_theta[idx])
eq_y = eq_radius[idx] * np.cos(eq_theta[idx])

#nrow = len(model_names)
nrow = 2
model_names = ['beta', 'alpha']
subplot_height = height/nrow

for irow in range(nrow):
  ax = fig.add_axes([0.1,0.1+irow*subplot_height,width,subplot_height], aspect='equal')
  ax.axis('off')
  
  # colorbar axis
  cax = fig.add_axes([0.12+width, 0.1+irow*subplot_height, 0.01, 0.8*subplot_height])

  # contourfill dlnvs
  tag = model_names[irow]
  cmap = plt.cm.get_cmap("jet_r")
  zz = model[tag]*100

  #levels = np.concatenate((np.arange(-6,0,1), np.arange(1,6.1,1)))
  #cs = ax.contour(xx, yy, zz, levels=levels, colors=('k',), linewidths=(0.1,))
  #plt.clabel(cs, fmt='%2.1f', colors='k', fontsize=5)

  levels = np.concatenate((np.arange(-6,0,0.5), np.arange(0.5,6.1,0.5)))
  cs = ax.contourf(xx, yy, zz, cmap=cmap, levels=levels, extend="both")
  cs.cmap.set_over('black')
  cs.cmap.set_under('purple')

  levels = np.arange(-6,6.1,1)
  cb = plt.colorbar(cs, cax=cax, ticks=levels, orientation="vertical")
  #cb.set_label('%', fontsize=10)
  cb.ax.set_title('(%)', fontsize=10)

  # overlay Vp/Vs ratio every dtheta degree
  dtheta = theta[1] - theta[0]
  ndtheta = int(np.round(np.deg2rad(1)/dtheta))

  for itheta in range(0, len(theta), ndtheta):
    theta1 = theta[itheta]
    kappa1 = kappa[:,itheta]

    dkappa1 = kappa1-1.8
    dkappa1 = dkappa1/np.max(np.abs(dkappa1))

    x1 = np.sin(theta1) * radius
    y1 = np.cos(theta1) * radius
    dz = radius[-1]*(ndtheta*dtheta)
    dx = np.cos(theta1)*dz
    dy = -1.0*np.sin(theta1)*dz
    ax.plot(x1+dx*dkappa1, y1+dy*dkappa1, 'r', lw=1.5)
    ax.plot(x1, y1, 'w', lw=1.0)
 
  # plot seismicity
  ax.plot(eq_x, eq_y, 'w+', markersize=2)

  ## colorbar for contourfill
  #if irow == 0:
  #  cax = fig.add_axes([0.3, 0.07, 0.4, 0.01])
  #  cb = plt.colorbar(cs, cax=cax, orientation="horizontal")
  #  cb.ax.tick_params(labelsize=8)
  #  cb.set_label('% REF', fontsize=10)
  #  #l, b, w, h = ax.get_position().bounds
  #  #ll, bb, ww, hh = CB.ax.get_position().bounds
  #  #cb.ax.set_position([ll, b + 0.1*h, ww, h*0.8])
 
  # mark certain depths 
  #for depth in [40, 220, 410, 670, 900]:
  for depth in [40, 220, 410]:
    x = np.sin(theta) * (R_earth_km - depth)
    y = np.cos(theta) * (R_earth_km - depth)
    ax.plot(x, y, 'k', lw=0.5)
    ax.text(x[-1],y[-1], str(depth), horizontalalignment='left', verticalalignment='top', fontsize=8)

  # surface marker
  marker_theta = np.linspace(theta[0], theta[-1], nmarker)
  x = np.sin(marker_theta) * radius[-1]
  y = np.cos(marker_theta) * radius[-1]
  ax.plot(x, y, 'ko', markersize=4, clip_on=False)
  x1 = np.sin(theta[-1]) * radius[-1]
  x2 = np.sin(theta[-2]) * radius[-1]
  y1 = np.cos(theta[-1]) * radius[-1]
  y2 = np.cos(theta[-2]) * radius[-1]
  xlim = ax.get_xlim()
  length = xlim[1]-xlim[0]
  ax.arrow(x1, y1, x1-x2, y1-y2, head_width=0.02*length, head_length=0.02*length, fc='k', ec='k', clip_on=False)

  # text model tag
  xlim = ax.get_xlim()
  ylim = ax.get_ylim()
  ax.text(xlim[0], ylim[1], tag, 
      horizontalalignment='left',
      verticalalignment='top', 
      fontsize=14)
  
#------ save figure
#plt.show()
plt.savefig(out_fig, format='pdf')