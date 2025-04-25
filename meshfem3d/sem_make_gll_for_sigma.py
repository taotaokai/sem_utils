#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" smoothing specfem3D GLL model in the lateral and vertical directions using Gaussian window of spatially varying size 

ECEF x,y,z coordinates

"""
import sys
import time

import numpy as np
from scipy.io import FortranFile

import pyproj

from meshfem3d_utils import sem_mesh_read
from meshfem3d_constants import NGLLX,NGLLY,NGLLZ,R_EARTH

#====== parameters
nproc = int(sys.argv[1])
mesh_dir = str(sys.argv[2]) # <mesh_dir>/proc******_external_mesh.bin
out_dir = str(sys.argv[3])

GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

lat_min = 22
lat_max = 23
#sigmaH_lat_min = 30 # km
sigmaH_lat_min = 5 # km
sigmaH_lat_max = 5 # km

depth_min = 150
depth_max = 190
sigmaV_depth_min = 1 #km
sigmaV_depth_max = 0 #km

#====== smooth each target mesh slice
for iproc in range(nproc):

  print("====== ", iproc)

  #--- read in SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir, iproc)
  mesh = sem_mesh_read(mesh_file)

  nspec = mesh['nspec']
  gll_dims = mesh['gll_dims']

  ibool = mesh['ibool']
  xyz_glob = mesh['xyz_glob'] * R_EARTH

  lon_glob, lat_glob, height_glob = pyproj.transform(ecef, lla, xyz_glob[0,:], xyz_glob[1,:], xyz_glob[2,:])

  depth_glob = -1.0 * height_glob / 1000.0 # km

  # sigmaH
  sigmaH = np.zeros(xyz_glob.shape[1]) 

  idx = lat_glob < lat_min
  sigmaH[idx] = sigmaH_lat_min

  idx = lat_glob >= lat_max
  sigmaH[idx] = sigmaH_lat_max

  idx = (lat_glob >= lat_min) & (lat_glob < lat_max)
  coef = (lat_glob[idx] - lat_min) / (lat_max - lat_min)
  sigmaH[idx] = (1.0 - coef)*sigmaH_lat_min + coef*sigmaH_lat_max
 
  # sigmaV
  sigmaV = np.zeros(xyz_glob.shape[1]) 

  idx = depth_glob < depth_min
  sigmaV[idx] = sigmaV_depth_min

  idx = depth_glob >= depth_max
  sigmaV[idx] = sigmaV_depth_max
  sigmaH[idx] = 0.0 # also make sigmaH zero below depth_max

  idx = (depth_glob >= depth_min) & (depth_glob < depth_max)
  coef = (depth_glob[idx] - depth_min) / (depth_max - depth_min)
  sigmaV[idx] = (1.0 - coef)*sigmaV_depth_min + coef*sigmaV_depth_max
  sigmaH[idx] = (1.0 - coef)*sigmaH[idx]

  #--- output
  out_file = "%s/proc%06d_reg1_sigmaH.bin"%(out_dir, iproc)
  with FortranFile(out_file, 'w') as f:
    f.write_record(np.array(sigmaH[ibool.ravel(order='F')-1], dtype='f4'))

  out_file = "%s/proc%06d_reg1_sigmaV.bin"%(out_dir, iproc)
  with FortranFile(out_file, 'w') as f:
    f.write_record(np.array(sigmaV[ibool.ravel(order='F')-1], dtype='f4'))
