#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" smoothing specfem3D GLL model in the lateral and vertical directions using Gaussian window of spatially varying size 

ECEF x,y,z coordinates

"""
import sys
import time

import numpy as np
from scipy import spatial
from scipy.io import FortranFile

import random

from mpi4py import MPI

from meshfem3d_utils import sem_mesh_read, sem_mesh_get_vol_gll
from meshfem3d_constants import *

from smooth_gauss_cap import smooth_gauss_cap 

#====== parameters
nproc = int(sys.argv[1])
mesh_dir = str(sys.argv[2]) # <mesh_dir>/proc******_external_mesh.bin
model_dir = str(sys.argv[3]) # <model_dir>/proc******_<model_name>.bin
model_names = str(sys.argv[4]) # comma delimited e.g. vp,vs,rho,qmu,qkappa

sigma_dir = str(sys.argv[5]) # <sigma_dir>/proc******_<sigmaH/R_tag>.bin
sigmaH_tag = str(sys.argv[6]) # tag for GLL files of sigmaH value (horizontal smoothing length in KM)
sigmaV_tag = str(sys.argv[7]) # tag for GLL files of sigmaV value (vertical/radial smoothing length in KM) 

out_dir = str(sys.argv[8])

#--- merge regions
#idoubling_merge = []
#In SETibet case, since I use a velocity gradient across Moho and no mesh boundary at Moho, treat IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80 as the same region
#idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220]
idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80]

#--- model names
model_names = model_names.split(',')
nmodel = len(model_names)

#====== smooth each target mesh slice
comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

if mpi_size < nproc:
  raise Exception("mpi_size must larger than nproc!")

if mpi_rank == 0:
  tic = time.time()

#--- read in SEM mesh slice
if mpi_rank < nproc:
  iproc_slice = mpi_rank
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir, iproc_slice)
  mesh_slice = sem_mesh_read(mesh_file)
  
  nspec_slice = mesh_slice['nspec']
  nglob_slice = mesh_slice['nglob']
  gll_dims_slice = mesh_slice['gll_dims']
  idoubling_slice = mesh_slice['idoubling']
  ibool_slice = mesh_slice['ibool']
  
  # merge regions if required
  idx_merge = np.zeros(nspec_slice, dtype='bool')
  for ii in idoubling_merge:
    idx_merge = idx_merge | (idoubling_slice == ii)
  idoubling_slice[idx_merge] = IFLAG_DUMMY
  
  xyz_elem_slice = mesh_slice['xyz_elem']
  xyz_glob_slice = mesh_slice['xyz_glob']
  
  ##DEBUG again the problem of non-contiguous array
  #if mpi_rank == 0:
  #  print(xyz_elem_slice.flags['C_CONTIGUOUS']) #FALSE
  #  print(xyz_glob_slice.flags['C_CONTIGUOUS']) #TRUE
  
  #--- read in smoothing length
  sigmaH_file = "%s/proc%06d_reg1_%s.bin"%(sigma_dir, iproc_slice, sigmaH_tag)
  with FortranFile(sigmaH_file, 'r') as f:
    # note: must use fortran convention when reshape to N-D array!!!
    sigmaH_gll_slice = np.ascontiguousarray(np.reshape(f.read_ints(dtype='f4'), gll_dims_slice, order='F'))
  
  sigmaV_file = "%s/proc%06d_reg1_%s.bin"%(sigma_dir, iproc_slice, sigmaV_tag)
  with FortranFile(sigmaV_file, 'r') as f:
    # note: must use fortran convention when reshape to N-D array!!!
    sigmaV_gll_slice = np.ascontiguousarray(np.reshape(f.read_ints(dtype='f4'), gll_dims_slice, order='F'))
  
  # non-dimensionalize as SEM
  sigmaH_gll_slice /= R_EARTH_KM
  sigmaV_gll_slice /= R_EARTH_KM
  
  #--- determine element size (approximately)
  max_gll_dist_slice = np.zeros(nspec_slice)
  min_gll_dist_slice = np.zeros(nspec_slice)
  for ispec in range(nspec_slice):
    # distance between gll points and the central gll point 
    iglob1 = ibool_slice[:,:,:,ispec].ravel() - 1
    xyz_gll = xyz_glob_slice[:,iglob1]
    dist2 = np.sum((xyz_gll[:,:,None] - xyz_gll[:,None,:])**2, axis=0)
    np.fill_diagonal(dist2, np.nan)
    #dist = np.sum((xyz_elem_slice[:,ispec:ispec+1] - xyz_glob_slice[:,iglob1])**2, axis=0)**0.5
    max_gll_dist_slice[ispec] = np.nanmax(dist2)**0.5
    min_gll_dist_slice[ispec] = np.nanmin(dist2)**0.5
    # replace small sigma values with a fraction of minimum distance between gll points
    # to eliminate the issue from diving zero in smooth_gauss_cap.f90
    idx = sigmaH_gll_slice[:,:,:,ispec] < 1.0e-4*min_gll_dist_slice[ispec]
    sigmaH_gll_slice[idx,ispec] = 1.0e-4*min_gll_dist_slice[ispec]
    idx = sigmaV_gll_slice[:,:,:,ispec] < 1.0e-4*min_gll_dist_slice[ispec]
    sigmaV_gll_slice[idx,ispec] = 1.0e-4*min_gll_dist_slice[ispec]
  
  # search neighbouring contrib elements
  search_radius2_slice = (5*np.maximum(np.amax(sigmaH_gll_slice, axis=(0,1,2)), np.amax(sigmaV_gll_slice, axis=(0,1,2))) + max_gll_dist_slice)**2
  search_radius2_slice = np.ascontiguousarray(np.array(search_radius2_slice, dtype=np.float32))
  #if mpi_rank == 0:
  #  print(search_radius2_slice.dtype)
  #comm.Barrier()
  #sys.exit(-1)
  
  #tree_elem_slice = spatial.cKDTree(np.column_stack((xyz_elem_slice[0,:],xyz_elem_slice[1,:],xyz_elem_slice[2,:])))
  
  #---- get gll integration weights
  vol_gll_slice = sem_mesh_get_vol_gll(mesh_slice)
  
  #---- read model values of the sliceuting mesh slice
  model_gll_slice = np.empty((nmodel,)+gll_dims_slice, dtype=np.float32) # use float32 and bcast as MPI.FLOAT
  for imodel in range(nmodel):
    model_tag = model_names[imodel]
    model_file = "%s/proc%06d_reg1_%s.bin"%(model_dir, iproc_slice, model_tag)
    with FortranFile(model_file, 'r') as f:
      # note: must use fortran convention when reshape to N-D array!!! 
      model_gll_slice[imodel,:,:,:,:] = np.reshape(f.read_ints(dtype='f4'), gll_dims_slice, order='F')

comm.Barrier()
if mpi_rank == 0:
  print("====== finish reading mesh/model ", time.time()-tic)
  sys.stdout.flush()

#====== loop over each target mesh slice
#for iproc_target in range(nproc):
for iproc_target in range(266,nproc):
#for iproc_target in range(240,nproc):
#for iproc_target in range(346,nproc):
#for iproc_target in range(246):

  if mpi_rank == 0:
    print("====== target ", iproc_target, time.time()-tic)
    sys.stdout.flush()

  #--- broadcast target mesh slice
  nspec_target = None
  nglob_target = None
  gll_dims_target = None

  if mpi_rank == iproc_target:
    nspec_target = nspec_slice
    nglob_target = nglob_slice
    gll_dims_target = gll_dims_slice

  nspec_target = comm.bcast(nspec_target, root=iproc_target)
  nglob_target = comm.bcast(nglob_target, root=iproc_target)
  gll_dims_target = comm.bcast(gll_dims_target, root=iproc_target)

  ibool_target = np.empty(gll_dims_target, dtype=np.int32)
  idoubling_target = np.empty(nspec_target, dtype=np.int32)
  xyz_glob_target = np.empty((3,nglob_target), dtype=np.float32)
  xyz_elem_target = np.empty((3,nspec_target), dtype=np.float32)
  search_radius2_target = np.empty(nspec_target, dtype=np.float32)
  sigmaH_gll_target = np.empty(gll_dims_target, dtype=np.float32)
  sigmaV_gll_target = np.empty(gll_dims_target, dtype=np.float32)

  if mpi_rank == iproc_target:
    ibool_target[:] = ibool_slice
    idoubling_target[:] = idoubling_slice
    xyz_glob_target[:] = xyz_glob_slice
    xyz_elem_target[:] = xyz_elem_slice
    search_radius2_target[:] = search_radius2_slice
    sigmaH_gll_target[:] = sigmaH_gll_slice
    sigmaV_gll_target[:] = sigmaV_gll_slice

    # only needed on iproc_target
    weight_model_gll_target = np.zeros((nmodel,)+gll_dims_target)
    weight_gll_target = np.zeros(gll_dims_target)

  comm.Bcast(ibool_target, root=iproc_target)
  comm.Bcast(idoubling_target, root=iproc_target)
  comm.Bcast(xyz_glob_target, root=iproc_target)
  comm.Bcast(xyz_elem_target, root=iproc_target)
  comm.Bcast(search_radius2_target, root=iproc_target)
  comm.Bcast(sigmaH_gll_target, root=iproc_target)
  comm.Bcast(sigmaV_gll_target, root=iproc_target)

  if mpi_rank == 0:
    print("finish bcast target mesh ", time.time()-tic)
    sys.stdout.flush()

  #====== loop over each contribution mesh slice
  for iproc_contrib in range(nproc):

    if mpi_rank == 0:
      print("--- contrib ", iproc_contrib)
      sys.stdout.flush()

    #--- broadcast contrib mesh slice
    nspec_contrib = None
    nglob_contrib = None
    gll_dims_contrib = None

    if mpi_rank == iproc_contrib:
      nspec_contrib = nspec_slice
      nglob_contrib = nglob_slice
      gll_dims_contrib = gll_dims_slice

    nspec_contrib = comm.bcast(nspec_contrib, root=iproc_contrib)
    nglob_contrib = comm.bcast(nglob_contrib, root=iproc_contrib)
    gll_dims_contrib = comm.bcast(gll_dims_contrib, root=iproc_contrib)

    idoubling_contrib = np.empty(nspec_contrib, dtype=np.int32)
    xyz_elem_contrib = np.empty((3,nspec_contrib), dtype=np.float32)

    if mpi_rank == iproc_contrib:
      idoubling_contrib[:] = idoubling_slice
      xyz_elem_contrib[:] = xyz_elem_slice

    comm.Bcast(idoubling_contrib, root=iproc_contrib)
    comm.Bcast(xyz_elem_contrib, root=iproc_contrib)

    if mpi_rank == 0:
      print("finish bcast contrib mesh ", time.time()-tic)
      sys.stdout.flush()

    #--- select contrib elements based on distance and idoubling to narrow down the elements that need explicit calculation

    # index contrib elements within a distance threshold
    dist2_elem_target_contrib_local = np.sum((xyz_elem_target[:,mpi_rank::mpi_size,None] - xyz_elem_contrib[:,None,:])**2, axis=0)
    idx_dist_local = dist2_elem_target_contrib_local < search_radius2_target[mpi_rank::mpi_size,None]

    # only smooth over contrib elements of the same idoubling
    idx_idoubling_local = (idoubling_target[mpi_rank::mpi_size,None] - idoubling_contrib[None,:]) == 0

    # final index of selected contrib elements
    idx_select_local = idx_dist_local & idx_idoubling_local

    # check if no contrib elements found
    nlocal = np.sum(idx_select_local)
    nlocal_sum = comm.allreduce(nlocal, op=MPI.SUM)
    if mpi_rank == 0:
      print("number of contrib elements found ", nlocal_sum, time.time()-tic)
      sys.stdout.flush()
    if nlocal_sum == 0:
      continue

    #--- bcast affected list of ispec_target and its associate ispec_contrib
    ispec_target_local = np.arange(mpi_rank,nspec_target,mpi_size, dtype=np.int32)
    nlocal = len(ispec_target_local)
    ispec_target_list_local = np.broadcast_to(ispec_target_local[:,None], (nlocal,nspec_contrib))[idx_select_local]
    ispec_contrib_list_local = np.broadcast_to(np.arange(nspec_contrib, dtype=np.int32), (nlocal,nspec_contrib))[idx_select_local]

    nlocal = ispec_target_list_local.size
    counts = np.empty(mpi_size, dtype=np.int32)
    comm.Allgather(np.array(nlocal,dtype=np.int32), counts)
    displs = np.cumsum(np.r_[([0,], counts[0:-1])])
    
    count_all = sum(counts)
    ispec_target_list = np.empty(count_all, dtype=np.int32)
    ispec_contrib_list = np.empty(count_all, dtype=np.int32)
    comm.Allgatherv([ispec_target_list_local, counts[mpi_rank]], [ispec_target_list, counts, displs, MPI.INT32_T])
    comm.Allgatherv([ispec_contrib_list_local, counts[mpi_rank]], [ispec_contrib_list, counts, displs, MPI.INT32_T])

    #comm.Gatherv([ispec_target_list_local, counts[mpi_rank]], [ispec_target_list, counts, displs, MPI.INT32_T])
    #comm.Gatherv([ispec_contrib_list_local, counts[mpi_rank]], [ispec_contrib_list, counts, displs, MPI.INT32_T])
    #if mpi_rank == iproc_target:
    #  ind = np.random.permutation(count_all)
    #  ispec_target_list = ispec_target_list[ind]
    #  ispec_contrib_list = ispec_contrib_list[ind]
    #comm.Bcast(ispec_target_list, root=iproc_target)
    #comm.Bcast(ispec_contrib_list, root=iproc_target)

    #--- bcast model
    ibool_contrib = np.empty(gll_dims_contrib, dtype=np.int32)
    xyz_glob_contrib = np.empty((3,nglob_contrib), dtype=np.float32)
    model_gll_contrib = np.empty((nmodel,)+gll_dims_contrib, dtype=np.float32)
    vol_gll_contrib = np.empty(gll_dims_contrib, dtype=np.float64)

    if mpi_rank == iproc_contrib:
      ibool_contrib[:] = ibool_slice
      xyz_glob_contrib[:] = xyz_glob_slice
      model_gll_contrib[:] = model_gll_slice
      vol_gll_contrib[:] = vol_gll_slice

    comm.Bcast(ibool_contrib, root=iproc_contrib)
    comm.Bcast(xyz_glob_contrib, root=iproc_contrib)
    comm.Bcast(model_gll_contrib, root=iproc_contrib)
    comm.Bcast(vol_gll_contrib, root=iproc_contrib)

    if mpi_rank == 0:
      print("finish bcast contrib model ", time.time()-tic)
      sys.stdout.flush()

    #--- gather contributions for each target gll point
    weight_model_gll_local = np.zeros((nmodel,)+gll_dims_target)
    weight_gll_local = np.zeros(gll_dims_target)

    ##NOTE if throw large amount of data to smooth_gauss_cap, it will be extremely slow. Contrary to what I expect.
    #ispec_target_local = ispec_target_list[mpi_rank::mpi_size]
    #ispec_contrib_local = ispec_contrib_list[mpi_rank::mpi_size]
    #nlocal = len(ispec_target_local)
    #iglob_target = ibool_target[:,:,:,ispec_target_local].ravel() - 1
    #xyz_gll_target = xyz_glob_target[:,iglob_target]
    #sigmaH = sigmaH_gll_target[:,:,:,ispec_target_local].ravel()
    #sigmaV = sigmaV_gll_target[:,:,:,ispec_target_local].ravel()
    #iglob_contrib = ibool_contrib[:,:,:,ispec_contrib_local].ravel() - 1
    #xyz_gll_contrib = xyz_glob_contrib[:,iglob_contrib]
    #model_gll = model_gll_contrib[:,:,:,:,ispec_contrib_local].reshape((nmodel,-1))
    #vol_gll = vol_gll_contrib[:,:,:,ispec_contrib_local].ravel()
    #weight_model_val, weight_val = smooth_gauss_cap(xyz_gll_target, xyz_gll_contrib, vol_gll, model_gll, sigmaH, sigmaV) 
    #weight_model_gll_local[:,:,:,:,ispec_target_local] += weight_model_val.reshape((nmodel,NGLLX,NGLLY,NGLLZ,nlocal))
    #weight_gll_local[:,:,:,ispec_target_local] += weight_val.reshape((NGLLX,NGLLY,NGLLZ,nlocal))

    ##NOTE I do this try to throw more data into smooth_gauss_cap, but really does not help with the speed.
    #ispec_target_unique, ind = np.unique(ispec_target_local, return_inverse=True)
    #for icount in range(len(ispec_target_unique)):
      #ispec_target = ispec_target_unique[icount]
      #ispec_contrib = ispec_contrib_local[ind==icount]

    for icount in range(mpi_rank,count_all,mpi_size):

      ispec_target = ispec_target_list[icount]
      ispec_contrib = ispec_contrib_list[icount]

      iglob_target = ibool_target[:,:,:,ispec_target].ravel() - 1
      xyz_gll_target = xyz_glob_target[:,iglob_target]
      sigmaH = sigmaH_gll_target[:,:,:,ispec_target].ravel()
      sigmaV = sigmaV_gll_target[:,:,:,ispec_target].ravel()

      iglob_contrib = ibool_contrib[:,:,:,ispec_contrib].ravel() - 1
      xyz_gll_contrib = xyz_glob_contrib[:,iglob_contrib]

      model_gll = model_gll_contrib[:,:,:,:,ispec_contrib].reshape((nmodel,-1))
      vol_gll = vol_gll_contrib[:,:,:,ispec_contrib].ravel()

      weight_model_val, weight_val = smooth_gauss_cap(xyz_gll_target, xyz_gll_contrib, vol_gll, model_gll, sigmaH, sigmaV) 

      weight_model_gll_local[:,:,:,:,ispec_target] += weight_model_val.reshape((nmodel,NGLLX,NGLLY,NGLLZ))
      weight_gll_local[:,:,:,ispec_target] += weight_val.reshape((NGLLX,NGLLY,NGLLZ))

    #print("# finish loop ilocal_target rank,len,elapse_time =  ", mpi_rank, len(ilocal_target_list), time.time()-tic)
    #sys.stdout.flush()
    #end for ispec_target in ispec_target_list[mpi_rank::mpi_size]:

    # gather results
    if mpi_rank == iproc_target:
      weight_model_gll_local_sum = np.zeros((nmodel,)+gll_dims_target)
      weight_gll_local_sum = np.zeros(gll_dims_target)
    else:
      weight_model_gll_local_sum = None
      weight_gll_local_sum = None

    #print("rank, max(weight_gll_local) = ", mpi_rank, np.max(weight_gll_local))

    comm.Reduce(weight_model_gll_local, weight_model_gll_local_sum, op=MPI.SUM, root=iproc_target)
    comm.Reduce(weight_gll_local, weight_gll_local_sum, op=MPI.SUM, root=iproc_target)

    if mpi_rank == iproc_target:
      print("max(weight_gll_local_sum) = ", np.max(weight_gll_local_sum))
      print("max(weight_model_gll_local_sum) = ", np.max(weight_model_gll_local_sum))
      sys.stdout.flush()
      weight_model_gll_target += weight_model_gll_local_sum
      weight_gll_target += weight_gll_local_sum

    if mpi_rank == 0:
      print("finish contrib mesh ", time.time()-tic)
      sys.stdout.flush()

    comm.Barrier()

  #END====== loop over each contributing mesh slice

  if mpi_rank == iproc_target:
    #--- weighted average model values on each target gll point
    weight_model_gll_target /= weight_gll_target
 
    #--- output model gll file
    for imodel in range(nmodel):
      model_tag = model_names[imodel]
      out_file = "%s/proc%06d_reg1_%s.bin"%(out_dir, iproc_target, model_tag)
      with FortranFile(out_file, 'w') as f:
        out_data = np.ravel(weight_model_gll_target[imodel,:], order='F') # Fortran's column-major convention
        f.write_record(np.array(out_data, dtype='f4'))

  comm.Barrier()
#END loop over each target mesh slice    
