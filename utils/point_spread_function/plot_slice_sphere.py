#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Plot .nc files generated by xsem_slice_sphere
"""
import sys
import numpy as np

from netCDF4 import Dataset

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#====== utility functions
def round_to_1(x):
  """round a number to one significant figure
  """
  return np.round(x, -int(np.floor(np.log10(np.abs(x)))))

# round to n significant figure
round_to_n = lambda x, n: np.round(x, -int(np.floor(np.log10(x))) + (n - 1)) 

#------ user inputs
nc_file = sys.argv[1]
title = sys.argv[2]
out_fig = sys.argv[3]

model_names = ['beta_Hdm_precond', 'alpha_Hdm_precond', 'dlnv']

#-- plot controls
# map plot
map_lat_center = 38.5
map_lon_center = 118.0
map_lat_min= 10
map_lat_max= 55
map_lon_min= 90
map_lon_max= 160
map_parallels = np.arange(0.,81,10.)
map_meridians = np.arange(0.,351,10.)

# model colormap
# model plot region
width = 0.8
height = 0.9

#------ read nc files 
fh = Dataset(nc_file, mode='r')

lats = fh.variables['latitude'][:]
lons = fh.variables['longitude'][:]

slice_depth = fh.variables['depth'][:]

model = {}
for tag in model_names:
  model[tag] = fh.variables[tag][:]

#------ plot map and xsection surface trace and marker
fig = plt.figure(figsize=(8.5,11))

# add title
ax = fig.add_axes([0.5, 0.07+height, 0.1, 0.1])
ax.axis('off')
ax.text(0, 0, title, 
    horizontalalignment='center',
    verticalalignment='top',
    fontsize=16)

lons2, lats2 = np.meshgrid(lons, lats, indexing='ij')
nrow = len(model_names)
subplot_height = height/nrow


zz_dlnv = np.transpose(model['dlnv']) * -100.0

for irow in range(nrow):
  origin_x = 0.1
  origin_y = 0.05+irow*subplot_height
  ax = fig.add_axes([origin_x, origin_y, width, 0.8*subplot_height])
  
  model_tag = model_names[irow]

  m = Basemap(ax=ax, projection='tmerc', resolution='l',
      llcrnrlat=map_lat_min, llcrnrlon=map_lon_min, 
      urcrnrlat=map_lat_max, urcrnrlon=map_lon_max,
      lat_0=map_lat_center, lon_0=map_lon_center)
  m.drawcoastlines(linewidth=0.2)
  m.drawcountries(linewidth=0.2)
  m.drawparallels(map_parallels, linewidth=0.1, labels=[1,0,0,0], fontsize=8)
  m.drawmeridians(map_meridians, linewidth=0.1, labels=[0,0,0,1], fontsize=8)
  
  # contourf model
  xx, yy = m(lons2, lats2)
  zz = np.transpose(model[model_tag])

  if model_tag == 'dlnv':
    #z_max = round_to_1(np.max(np.abs(zz)))
    #dz = 2.0*z_max/10
    #levels = np.arange(-z_max, z_max+dz/2, dz)
    zz = zz * -100.0 # use percentage

    #levels = np.concatenate((np.arange(-6,0,1), np.arange(1,6.1,1)))
    #cs = m.contour(xx, yy, zz, levels=levels, colors=('k',), linewidths=(0.1,))
    #plt.clabel(cs, fmt='%1.0f', colors='k', fontsize=3)

    levels = np.linspace(-1,1,100)
    cmap = plt.cm.get_cmap("jet_r")
    cs = m.contourf(xx, yy, zz, cmap=cmap, levels=levels, extend="both")
    cs.cmap.set_over('black')
    cs.cmap.set_under('purple')
    # colorbar for contourfill
    cb = m.colorbar(cs,location='right',ticks=np.arange(-1,1.1,0.5), pad="5%")
    cb.ax.set_title('(%)', fontsize=10)
    #cb.set_label('%')
    ax.set_title("{:s}".format(model_tag))

  elif model_tag in ['alpha_Hdm_precond', 'beta_Hdm_precond']:
    #z_max = round_to_1(np.max(np.abs(zz)))
    #dz = 2.0*z_max/10
    #levels = np.arange(-z_max, z_max+dz/2, dz)
    #norm_amp = np.max(zz)
    #norm_amp = 4.0e-6
    norm_amp = 10
    zz /= norm_amp

    levels = [-np.exp(-0.5), np.exp(-0.5),]
    cs = m.contour(xx, yy, zz_dlnv, levels=levels, colors=('w',), linewidths=(0.1,))
    #plt.clabel(cs, fmt='%1.0f', colors='k', fontsize=3)

    levels = np.linspace(-1,1,100)
    cmap = plt.cm.get_cmap("jet_r")
    cs = m.contourf(xx, yy, zz, cmap=cmap, levels=levels, extend="both")
    cs.cmap.set_over('black')
    cs.cmap.set_under('purple')
    # colorbar for contourfill
    cb = m.colorbar(cs,location='right',ticks=np.arange(-1,1.1,0.5), pad="5%")
    cb.ax.set_title('x %7.1E'%(norm_amp), fontsize=10)
    #cb.set_label('%')
    ax.set_title("{:s}".format(model_tag))

  #elif model_tag in ['kappa']:
  #  dz = (np.max(zz) - np.min(zz))/10
  #  dz = round_to_1(dz)
  #  levels = np.arange(np.min(zz), np.max(zz)+dz/2, dz)
  #  cmap = plt.cm.get_cmap("jet")
  #  cs = m.contourf(xx, yy, zz, cmap=cmap, levels=levels, extend="both")
  #  cs.cmap.set_over('purple')
  #  cs.cmap.set_under('black')
  #  # colorbar for contourfill
  #  cb = m.colorbar(cs,location='right',pad="5%", format="%.2f")
  #  #cb.set_label('% mean')
  #  ax.set_title("{:s}".format(model_tag))

  #elif model_tag in ['gamma', 'eps']:
  #  dz = (np.max(zz) - np.min(zz))/10
  #  if dz < np.finfo(np.float).eps:
  #    levels = (0, 1.0)
  #  else:
  #    z_max = round_to_1(np.max(np.abs(zz)))
  #    dz = 2.0*z_max/10
  #    levels = np.arange(-z_max, z_max+dz/2, dz)
  #    #dz = round_to_1(dz)
  #    #levels = np.arange(np.min(zz), np.max(zz)+dz/2, dz)
  #  cmap = plt.cm.get_cmap("jet")
  #  cs = m.contourf(xx, yy, zz, cmap=cmap, levels=levels, extend="both")
  #  cs.cmap.set_over('purple')
  #  cs.cmap.set_under('black')
  #  # colorbar for contourfill
  #  cb = m.colorbar(cs,location='right',pad="5%", format="%.2f")
  #  #cb.set_label('% mean')
  #  ax.set_title("{:s}".format(model_tag))
  else:
    raise Exception("unrecognized model {:s}".format(model_tag))

  #-- plot fault lines
  #fault_line_file = 'fault_lines.txt'
  #fault_lines = []
  #with open(fault_line_file, 'r') as f:
  #  lon = []
  #  lat = []
  #  for l in f.readlines():
  #    if not l.startswith('>'):
  #      x = l.split()
  #      lon.append(float(x[0]))
  #      lat.append(float(x[1]))
  #    else:
  #      fault_lines.append([lon, lat])
  #      lon = []
  #      lat = []
  #for l in fault_lines:
  #  x, y = m(l[0], l[1])
  #  ax.plot(x, y, 'k-', lw=0.05)

#  #--- plot seismicity
#  catalog_file = 'isc_d50km.txt'
#  with open(catalog_file, 'r') as f:
#    lines = [ l.split('|') for l in f.readlines() if not l.startswith('#') ] 
#    eq_lats = np.array([float(x[2]) for x in lines])
#    eq_lons = np.array([float(x[3]) for x in lines])
#    eq_deps = np.array([float(x[4]) for x in lines]) # km
#
#  eq_indx = np.abs(eq_deps - slice_depth) <= 10.0
#  x, y = m(eq_lons[eq_indx], eq_lats[eq_indx])
#  if slice_depth < 200:
#    markersize = 0.5
#  elif slice_depth < 500:
#    markersize = 1
#  else:
#    markersize = 2
#  ax.plot(x, y, 'w.', markersize=markersize)

#------ save figure
#plt.show()
plt.savefig(out_fig, format='pdf')