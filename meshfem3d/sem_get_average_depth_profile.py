#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys

import numpy as np
from scipy.io import FortranFile

import pyproj

from meshfem3d_utils import sem_mesh_read, sem_mesh_get_vol_gll, R_EARTH

#====== 
nproc = int(sys.argv[1])
mesh_dir = str(sys.argv[2]) # <mesh_dir>/proc******_external_mesh.bin
model_dir = str(sys.argv[3])
model_name = str(sys.argv[4])
out_dir = str(sys.argv[5])

depth_interp = np.hstack((np.arange(0,60,3), np.arange(60,101,5)))
half_depth_win = np.hstack((depth_interp[1]-depth_interp[0], np.diff(depth_interp)))/2.0
ndepth = depth_interp.size

GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

model_interp = np.zeros((nproc,ndepth))

for iproc in range(nproc):
  
  print("#proc = ", iproc)

  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir, iproc)
  mesh_data = sem_mesh_read(mesh_file)
  nspec = mesh_data['nspec']
  ibool = mesh_data['ibool']
  xyz_glob = mesh_data['xyz_glob']
  gll_dims = mesh_data['gll_dims']

  vol_gll = sem_mesh_get_vol_gll(mesh_data)

  lon, lat, height = pyproj.transform(ecef, lla, 
      xyz_glob[0,:]*R_EARTH, xyz_glob[1,:]*R_EARTH, xyz_glob[2,:]*R_EARTH)

  depth_gll = -1*height[ibool-1]/1000.0

  model_file = "%s/proc%06d_reg1_%s.bin"%(model_dir, iproc, model_name)
  with FortranFile(model_file, 'r') as f:
    model_gll = np.reshape(f.read_ints(dtype='f4'), gll_dims, order='F')
 
  for idep in range(ndepth):
    dep = depth_interp[idep]
    min_depth = dep - half_depth_win[idep]
    max_depth = dep + half_depth_win[idep]
    idx = (depth_gll >= min_depth) & (depth_gll <= max_depth)
    weight = 1.0 - np.abs(depth_gll[idx] - dep)/half_depth_win[idep]
    model_interp[iproc,idep] = np.sum(weight*vol_gll[idx]*model_gll[idx])/np.sum(weight*vol_gll[idx])

  out_file  = '%s/proc%06d_reg1_%s.txt'%(out_dir,iproc,model_name)
  with open(out_file, 'w') as f:
    for idep in range(ndepth):
      f.write("%12.5e  %12.5e\n"%(depth_interp[idep], model_interp[iproc,idep]))

model_avg = np.mean(model_interp,axis=0)
out_file  = '%s/reg1_average_%s.txt'%(out_dir,model_name)
with open(out_file, 'w') as f:
  for idep in range(ndepth):
    f.write("%12.5e  %12.5e\n"%(depth_interp[idep], model_avg[idep]))

