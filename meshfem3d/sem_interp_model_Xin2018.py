#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Merge Xin (2018) model within FWEA18
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
model_input_dir = str(sys.argv[1]) # horizontal model slices lon_lat_model consisting of "depth vp ..."
model_names = str(sys.argv[2]) # comma delimited e.g. vp,vs,rho,...

# mesh files for which to get interpolated values
nproc_target = int(sys.argv[3])
mesh_dir_target = str(sys.argv[4]) # <mesh_dir>/proc******_external_mesh.bin
model_dir_target = str(sys.argv[5]) # <model_dir>/proc******_<model_name>.bin

out_dir = str(sys.argv[6])

# model names
model_names = model_names.split(',')
nmodel = len(model_names)

# 3-D model grid of Xin2018 
lon_grid = np.hstack(([73,],np.arange(74,135.1,0.5), [136,]))
lat_grid = np.hstack(([17,],np.arange(18,53.1,0.5), [54,]))
depth_grid = np.array([0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 120, 150]) # km, relative to mean sea level

# merge depth (km), linear transition to FWEA18
merge_depth1 = 120
merge_depth2 = 150

# convert ECEF to Lon,Lat,Height
GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

#====== read in model file
nlon = lon_grid.size
nlat = lat_grid.size
ndepth = depth_grid.size
model_grid3 = np.empty((nlon,nlat,ndepth,nmodel))

for idep in range(ndepth):
  for imodel in range(nmodel):

    model_tag = model_names[imodel]
    model_file = "%s/Z_%s%d.txt"%(model_input_dir,model_tag[0:2],depth_grid[idep]) # only first two letters in model_tag, so vsv and vsh will be vs. 

    with open(model_file, 'r') as f:
      lines = [l.split() for l in f.readlines()]
    model_depth = np.array([float(l[2]) for l in lines])
    model_grid3[:,:,idep,imodel] = np.reshape(model_depth, (nlon,nlat))

#====== interpolate
for iproc_target in range(nproc_target):

  print("====== iproc_target ", iproc_target)
  sys.stdout.flush()

  # read in target SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_target, iproc_target)
  mesh_data_target = sem_mesh_read(mesh_file)
  nspec_target = mesh_data_target['nspec']
  ibool_target = mesh_data_target['ibool']
  xyz_glob_target = mesh_data_target['xyz_glob']

  # xyz points to interpolate
  # note: must use fortran convention when ravelling
  xyz_target = R_EARTH * xyz_glob_target[:,ibool_target.ravel(order='F')-1]

  # convert to lon,lat,ellipsoidal height
  llh_target = np.empty(xyz_target.shape)
  llh_target[0,:], llh_target[1,:], llh_target[2,:]  = pyproj.transform(ecef, lla, xyz_target[0,:], xyz_target[1,:], xyz_target[2,:])
  llh_target[2,:] = llh_target[2,:]/1000.0 # km
  llh_target = np.swapaxes(llh_target,0,1) # for interpn

  # convert to lon,lat, and depth below mean sea level
  lld_target = llh_target
  lld_target[:,2] = -1*(lld_target[:,2])

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

  # merge interpolation results into target model
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
