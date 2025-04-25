#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
from scipy.io import FortranFile
               
import pyproj

from meshfem3d_utils import sem_mesh_read, R_EARTH

#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
mesh_dir = str(sys.argv[3]) # <mesh_dir>/proc******_external_mesh.bin
vs_dir = str(sys.argv[4])
vs_tag = str(sys.argv[5])
vp_dir = str(sys.argv[6])
vp_tag = str(sys.argv[7])

GPS_ELLPS = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=GPS_ELLPS)
lla = pyproj.Proj(proj='latlong', ellps=GPS_ELLPS)

#====== read in gll file
for iproc in range(procnum_begin, procnum_end):

  print("#====== iproc", iproc)

  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir, iproc)
  mesh_data = sem_mesh_read(mesh_file)
  nspec = mesh_data['nspec']
  ibool = mesh_data['ibool']
  xyz_glob = mesh_data['xyz_glob']
  gll_dims = mesh_data['gll_dims']

  input_file = "%s/proc%06d_reg1_%s.bin"%(vs_dir, iproc, vs_tag)
  with FortranFile(input_file, 'r') as f:
    vs = np.reshape(f.read_ints(dtype='f4'), gll_dims, order='F')

  input_file = "%s/proc%06d_reg1_%s.bin"%(vp_dir, iproc, vp_tag)
  with FortranFile(input_file, 'r') as f:
    vp = np.reshape(f.read_ints(dtype='f4'), gll_dims, order='F')

  vp_vs_ratio = vp/vs

  ind = np.unravel_index(np.argmin(vp_vs_ratio, axis=None), vp_vs_ratio.shape)
  xyz = xyz_glob[:,ibool[ind]-1] * R_EARTH
  lon, lat, height = pyproj.transform(ecef, lla, xyz[0], xyz[1], xyz[2])
  print("min vp/vs, lon/lat/ellip_height = ", vp_vs_ratio[ind], lon, lat, height/1000.0)
