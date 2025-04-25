#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys

import numpy as np
from scipy.io import FortranFile

from meshfem3d_utils import sem_mesh_read, R_EARTH_KM

#====== 
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
mesh_dir = str(sys.argv[3]) # <mesh_dir>/proc******_external_mesh.bin
model_dir = str(sys.argv[4])
model_name = str(sys.argv[5])
DT = float(sys.argv[6])

for iproc in range(procnum_begin, procnum_end):

  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir, iproc)
  mesh_data = sem_mesh_read(mesh_file)
  nspec = mesh_data['nspec']
  ibool = mesh_data['ibool']
  xyz_glob = mesh_data['xyz_glob']
  gll_dims = mesh_data['gll_dims']
  
  model_file = "%s/proc%06d_reg1_%s.bin"%(model_dir, iproc, model_name)
  with FortranFile(model_file, 'r') as f:
    vel = np.reshape(f.read_ints(dtype='f4'), gll_dims, order='F')
  
  #--- determine element size (approximately)
  CFL = np.zeros(nspec)
  for ispec in range(nspec):
    iglob = ibool[:,:,:,ispec].ravel() - 1
    xyz_gll = xyz_glob[:,iglob]
    dist2 = np.sum((xyz_gll[:,:,None] - xyz_gll[:,None,:])**2, axis=0)**0.5*R_EARTH_KM
    np.fill_diagonal(dist2, np.nan)
    min_dist_gll = np.nanmin(dist2)
    CFL[ispec] = np.max(vel[:,:,:,ispec])*DT/min_dist_gll
  
  print("iproc, min/max vel = ", iproc, np.min(vel), np.max(vel))
  print("iproc, min/max CFL = ", iproc, np.min(CFL), np.max(CFL))
