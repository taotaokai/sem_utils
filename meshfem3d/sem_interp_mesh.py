#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys
import warnings
import time

import numpy as np
from scipy.io import FortranFile

from mpi4py import MPI

from meshfem3d_constants import NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA
from meshfem3d_constants import IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80,IFLAG_670_220,IFLAG_DUMMY
from gll_library import zwgljd, lagrange_poly

from meshfem3d_utils import sem_mesh_read, sem_locate_points_hex27

#====== parameters

# mesh files to interpolate from
nproc_source = int(sys.argv[1])
mesh_dir_source = str(sys.argv[2]) # <mesh_dir>/proc******_external_mesh.bin
model_dir_source = str(sys.argv[3]) # <model_dir>/proc******_<model_name>.bin

# mesh files for which to get interpolated values
nproc_target = int(sys.argv[4])
mesh_dir_target = str(sys.argv[5]) # <mesh_dir>/proc******_external_mesh.bin
#model_dir_target = str(sys.argv[6]) # <model_dir>/proc******_<model_name>.bin

model_names = str(sys.argv[6]) # comma delimited e.g. vp,vs,rho,qmu,qkappa
out_dir = str(sys.argv[7])

# merge regions
idoubling_merge = []
#In SETibet case, since I use a velocity gradient across Moho and no mesh boundary at Moho, treat IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80 as the same region
#idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220]
#idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80]

# model names
model_names = model_names.split(',')
nmodel = len(model_names)

#====== interpolate
comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

# GLL
xigll, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
yigll, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
zigll, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

#--- loop over each slice of target SEM mesh
for iproc_target in range(mpi_rank,nproc_target,mpi_size):

  print("====== iproc_target ", iproc_target)
  sys.stdout.flush()

  # read in target SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_target, iproc_target)
  mesh_data_target = sem_mesh_read(mesh_file)
  nspec_target = mesh_data_target['nspec']
  ibool_target = mesh_data_target['ibool']
  idoubling_target = mesh_data_target['idoubling']
  xyz_glob_target = mesh_data_target['xyz_glob']

  # merge regions if required
  idx_merge = np.zeros(nspec_target, dtype='bool')
  for ii in idoubling_merge:
    idx_merge = idx_merge | (idoubling_target == ii)
  idoubling_target[idx_merge] = IFLAG_DUMMY

  # xyz points to locate
  xyz_target = xyz_glob_target[:,ibool_target.ravel()-1]
  idoubling_ext = np.zeros(ibool_target.shape,dtype='int') + idoubling_target
  idoubling_ext = idoubling_ext.ravel()

  # array of final results
  npoints = xyz_target.shape[1]
  status_gll_target = np.zeros(npoints,dtype='int')
  status_gll_target[:] = -1
  misloc_gll_target = np.zeros(npoints)
  misloc_gll_target[:] = np.inf
  misratio_gll_target = np.zeros(npoints)
  model_gll_target = np.zeros((nmodel,npoints))

  # loop over each slice of source SEM mesh
  for iproc_source in range(nproc_source):

    print("iproc_target / iproc_source ", iproc_target, iproc_source)
    sys.stdout.flush()

    # read in source SEM mesh
    mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_source, iproc_source)
    mesh_data_source = sem_mesh_read(mesh_file)

    # merge regions if required
    idoubling_source = mesh_data_source['idoubling']
    idx_merge = np.zeros(mesh_data_source['nspec'], dtype='bool')
    for ii in idoubling_merge:
      idx_merge = idx_merge | (idoubling_source == ii)
    idoubling_source[idx_merge] = IFLAG_DUMMY

    # read in source model
    gll_dims = mesh_data_source['gll_dims']
    source_model_gll = np.zeros((nmodel,)+gll_dims)
    for imodel in range(nmodel):
      model_tag = model_names[imodel]
      model_file = "%s/proc%06d_reg1_%s.bin"%(model_dir_source, iproc_source, model_tag)
      with FortranFile(model_file, 'r') as f:
        # note: must use fortran convention when reshape to N-D array!!!
        source_model_gll[imodel,:,:,:,:] = np.reshape(f.read_ints(dtype='f4'), 
            gll_dims, order='F')

    # locate target points
    status_all, ispec_all, uvw_all, misloc_all, misratio_all = \
        sem_locate_points_hex27(mesh_data_source, xyz_target, idoubling_ext, kdtree_num_element=2.0, max_dist_ratio=2.0)

    # merge interpolation results of mesh slice (iproc_souce) into 
    # the final results based on misloc and status  

    # index selection for merge: (not located inside an element yet) and (located for the current mesh slice) and ( smaller misloc or located inside an element in this mesh slice )
    ii = (status_gll_target != 1) & (status_all != -1) & ( (misloc_all < misloc_gll_target) | (status_all == 1) )

    status_gll_target[ii] = status_all[ii]
    misloc_gll_target[ii] = misloc_all[ii]
    misratio_gll_target[ii] = misratio_all[ii]

    #hlagx = lagrange_poly(xigll, uvw_all[0,ii]).reshape((NGLLX,1,1,-1))
    #hlagy = lagrange_poly(yigll, uvw_all[1,ii]).reshape((1,NGLLY,1,-1))
    #hlagz = lagrange_poly(zigll, uvw_all[2,ii]).reshape((1,1,NGLLZ,-1))
    #model_gll_target[:,ii] = np.sum(source_model_gll[:,:,:,:,ispec_all[ii]]*hlagx*hlagy*hlagz, axis=(1,2,3))

    #FIXME replace with index slicing
    # status_gll_target[i] = status_all[ii]
    # ...
    # modify lagrange_poly(xigll,x) to handle x(n) 
    # ipoint_select = np.nonzero(ii)[0]
    #DONE

    #NOTE avoid too many for loops reduces computation time
    ##for ipoint in range(npoints):
    # slower than index slicing but use less memory 
    ipoint_select = np.nonzero(ii)[0]
    for ipoint in ipoint_select:
      #if (status_all[ipoint]==1 and status_gll_target[ipoint]==1):
      #  warnings.warn("point is located inside more than one element", 
      #      xyz_target[:,ipoint])
      # nothing to do if the point is already located inside an element
      # this means if multiple elements overlap (should not occur) we only take the first found element where the point locates inside
      #if status_gll_target[ipoint] == 1:
      #  continue
      #if (misloc_all[ipoint] > misloc_gll_target[ipoint]
      #    and status_all[ipoint]==1):
      #  warnings.warn("point located inside an element but with a larger misloc(loc/previous)", misloc_all['ipoint'], misloc_gll_target[ipoint])
      #if (misloc_all[ipoint] < misloc_gll_target[ipoint] 
      #    or status_all[ipoint]==1):
      #status_gll_target[ipoint] = status_all[ipoint]
      #misloc_gll_target[ipoint] = misloc_all[ipoint]
      #misratio_gll_target[ipoint] = misratio_all[ipoint]
      hlagx = lagrange_poly(xigll, uvw_all[0,ipoint])[:,0]
      hlagy = lagrange_poly(yigll, uvw_all[1,ipoint])[:,0]
      hlagz = lagrange_poly(zigll, uvw_all[2,ipoint])[:,0]
      model_gll_target[:,ipoint] = np.sum(source_model_gll[:,:,:,:,ispec_all[ipoint]]*hlagx[:,None,None]*hlagy[None,:,None]*hlagz[None,None,:], axis=(1,2,3))
        
  #end for loop over each slice of source SEM mesh

  # rehape results
  gll_dims = mesh_data_target['gll_dims']
  status_gll_target = np.reshape(status_gll_target, gll_dims)
  misloc_gll_target = np.reshape(misloc_gll_target, gll_dims)
  misratio_gll_target = np.reshape(misratio_gll_target, gll_dims)
  gll_dims = (nmodel,) + gll_dims
  model_gll_target = np.reshape(model_gll_target, gll_dims)

  # save interpolated model
  for imodel in range(nmodel):
    model_tag = model_names[imodel]
    model_file = "%s/proc%06d_reg1_%s.bin"%(out_dir, iproc_target, model_tag)
    with FortranFile(model_file, 'w') as f:
      f.write_record(np.array(np.ravel(model_gll_target[imodel,:,:,:,:], order='F'), dtype='f4'))

  # save misloc, status
  model_file = "%s/proc%06d_reg1_status.bin"%(out_dir, iproc_target)
  with FortranFile(model_file, 'w') as f:
    f.write_record(np.array(np.ravel(status_gll_target, order='F'), dtype='f4'))

  model_file = "%s/proc%06d_reg1_misloc.bin"%(out_dir, iproc_target)
  with FortranFile(model_file, 'w') as f:
    f.write_record(np.array(np.ravel(misloc_gll_target, order='F'), dtype='f4'))

# model_file = "%s/proc%06d_reg1_misratio.bin"%(out_dir, iproc_target)
# with FortranFile(model_file, 'w') as f:
#   f.write_record(np.array(np.ravel(misratio_gll_target, order='F'), dtype='f4'))
