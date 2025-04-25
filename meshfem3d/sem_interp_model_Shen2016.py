#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Merge Yang (2012) Tibet model within FWEA18
"""
import sys

import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import interp1d, interpn

from netCDF4 import Dataset
import pyproj

from meshfem3d_constants import R_EARTH
from meshfem3d_utils import sem_mesh_read

#====== parameters
model_input_dir = str(sys.argv[1]) # 1-D depth profiles lon_lat_model consisting of "depth vp ..."
model_names = str(sys.argv[2]) # comma delimited e.g. vp,vs,rho,...

# mesh files for which to get interpolated values
nproc_target = int(sys.argv[3])
mesh_dir_target = str(sys.argv[4]) # <mesh_dir>/proc******_external_mesh.bin
model_dir_target = str(sys.argv[5]) # <model_dir>/proc******_<model_name>.bin

out_dir = str(sys.argv[6])

# model names
model_names = model_names.split(',')
nmodel = len(model_names)

# topo
topoGRDfile = 'ETOPO1_Ice_g_smooth_b15km_I4m.grd'
GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

# 3-D grid of input 1-D profiles
lon_grid = np.arange(90,112,0.5)
lat_grid = np.arange(19,37,0.5)
#depth_grid = np.hstack((np.arange(-1,15,0.1), np.arange(15,149.8,0.2)))
depth_grid = np.arange(-1,190.1,0.5)

# merge depth (km), linear transition to FWEA18
merge_depth1 = 150
merge_depth2 = 190

#====== read in topo file
fh = Dataset(topoGRDfile,'r')

grd_lat = np.array(fh.variables['lat'][:])
grd_lon = np.array(fh.variables['lon'][:])
grd_topo = np.transpose(np.array(fh.variables['z'][:])) # (lon,lat)

#====== read in model file
nlon = lon_grid.size
nlat = lat_grid.size
ndepth = depth_grid.size
model_grid3 = np.empty((nlon,nlat,ndepth,nmodel))

print("====== Read models from ", model_input_dir)
for ilon in range(nlon):
  for ilat in range(nlat):

    model_file = "%s/%g_%g.mod"%(model_input_dir,lon_grid[ilon],lat_grid[ilat])

    print("--- lon lat model_file ", lon_grid[ilon],lat_grid[ilat], model_file)

    try:
      fh = open(model_file, 'r')
    except:
      fh = None
     
    lines = []
    if fh:
      lines = [l.split() for l in fh.readlines()]

    if len(lines) != 0:
      print("read model file ", model_file)
      model_depth = np.array([float(l[0]) for l in lines])
      ndepth = model_depth.size
      model_v = np.zeros((ndepth,nmodel))
      for imodel in range(nmodel):
        model_v[:,imodel] = np.array([float(l[1+imodel]) for l in lines])
      # add an additional negative depth point
      model_depth = np.hstack(([-1,], model_depth))
      model_v = np.vstack((model_v[0:1,:], model_v))
      f_interp = interp1d(model_depth, model_v, axis=0)
      model_grid3[ilon,ilat,:,:] = f_interp(depth_grid)
      fh.close()
    else:
      model_grid3[ilon,ilat,:,:] = np.nan

#====== interpolate
for iproc_target in range(nproc_target):
#for iproc_target in [134,]:

  print("====== iproc_target ", iproc_target)
  sys.stdout.flush()

  # read in target SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_target, iproc_target)
  mesh_data_target = sem_mesh_read(mesh_file)
  nspec_target = mesh_data_target['nspec']
  ibool_target = mesh_data_target['ibool']
  xyz_glob_target = mesh_data_target['xyz_glob']

  # xyz points to interpolate
  # note: must use fortran convention when ravelling to be consistent with model_interp/target
  xyz_target = R_EARTH * xyz_glob_target[:,ibool_target.ravel(order='F')-1]

  # convert to lon,lat,ellipsoidal height
  llh_target = np.empty(xyz_target.shape)
  llh_target[0,:], llh_target[1,:], llh_target[2,:]  = pyproj.transform(ecef, lla, xyz_target[0,:], xyz_target[1,:], xyz_target[2,:])
  llh_target[2,:] = llh_target[2,:]/1000.0 # km
  llh_target = np.swapaxes(llh_target,0,1) # for interpn
  #print(llh_target.shape)
  #print(llh_target[:,0:2].shape)
  topo_target = interpn((grd_lon,grd_lat), grd_topo, llh_target[:,0:2]) / 1000.0 #km

  # convert to lon,lat,depth below surface 
  lld_target = llh_target
  lld_target[:,2] = -1*(lld_target[:,2] - topo_target)
  #print(lld_target[0:10,:])

  # read in target model
  npoints = xyz_target.shape[1]
  model_target = np.zeros((npoints,nmodel))
  for imodel in range(nmodel):
    model_tag = model_names[imodel]
    gll_file = "%s/proc%06d_reg1_%s.bin"%(model_dir_target, iproc_target, model_tag)
    with FortranFile(gll_file, 'r') as f:
      #model_gll_target[imodel,:,:,:,:] = np.reshape(f.read_ints(dtype='f4'), gll_dims, order='F')
      model_target[:,imodel] = f.read_ints(dtype='f4')

  # interpolate
  model_interp = interpn((lon_grid,lat_grid,depth_grid), model_grid3, lld_target, bounds_error=False, fill_value=np.nan)
  #print(model_interp.shape)
  #print(model_interp[0:10])

  # merge interpolation results into target model

  #idx0 = ~(np.isnan(model_interp[:,0]))
  #model_target[idx0,:] = model_interp[idx0,:]

  print("len(isnan) = ", np.sum(np.isnan(model_interp[:,0])))

  idx0 = ~(np.isnan(model_interp[:,0]))
  idx1 = lld_target[:,2] < merge_depth1
  idx2 = (lld_target[:,2] >= merge_depth1) &  (lld_target[:,2] < merge_depth2)
  coef = (merge_depth2 - lld_target[:,2])/(merge_depth2 - merge_depth1)

  model_target[idx0&idx1,:] = model_interp[idx0&idx1,:]
  model_target[idx0&idx2,:] = model_interp[idx0&idx2,:]*coef[idx0&idx2,None] + (1 - coef[idx0&idx2,None])*model_target[idx0&idx2,:]

  # save interpolated model
  for imodel in range(nmodel):
    model_tag = model_names[imodel]
    gll_file = "%s/proc%06d_reg1_%s.bin"%(out_dir, iproc_target, model_tag)
    with FortranFile(gll_file, 'w') as f:
      f.write_record(np.array(model_target[:,imodel], dtype='f4'))
