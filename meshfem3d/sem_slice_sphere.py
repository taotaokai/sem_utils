#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys
import time

import numpy as np
from scipy.io import FortranFile
import pyproj
import xarray as xr
#from netCDF4 import Dataset
#from mpi4py import MPI

from meshfem3d_constants import NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA,R_EARTH
from gll_library import zwgljd, lagrange_poly
from meshfem3d_utils import sem_mesh_read, sem_locate_points_hex27

#====== parameters
nproc = int(sys.argv[1])
mesh_dir = str(sys.argv[2]) # <mesh_dir>/proc******_external_mesh.bin
model_dir = str(sys.argv[3]) # <model_dir>/proc******_<model_name>.bin
model_names = str(sys.argv[4]) # comma delimited e.g. vp,vs,rho,qmu,qkappa
model_units = str(sys.argv[5]) # comma delimited e.g. km.s-1,km.s-1,kg.m-3,count,count
min_lat = float(sys.argv[6])
max_lat = float(sys.argv[7])
nlat = int(sys.argv[8])
min_lon = float(sys.argv[9])
max_lon = float(sys.argv[10])
nlon = int(sys.argv[11])
depth_km = float(sys.argv[12]) # negative ellipsoidal height or ~ below mean sea level
out_file = str(sys.argv[13])

# model names
model_names = model_names.split(',')
model_units = model_units.split(',')
nmodel = len(model_names)

#--- create slice grid
grd_lon1 = np.linspace(min_lon, max_lon, nlon)
grd_lat1 = np.linspace(min_lat, max_lat, nlat)
grd_lon2, grd_lat2 = np.meshgrid(grd_lon1, grd_lat1, indexing='ij')
grd_alt2 = np.ones(grd_lon2.shape) * -1000.0 * depth_km
# convert (lon,lat,alt) to ECEF
ref_ellps = 'WGS84'
ecef = pyproj.Proj(proj='geocent', ellps=ref_ellps)
lla = pyproj.Proj(proj='latlong', ellps=ref_ellps)
xx, yy, zz = pyproj.transform(lla, ecef, grd_lon2, grd_lat2, grd_alt2)
xyz_interp = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())) / R_EARTH
npts_interp = xyz_interp.shape[1]

#====== interpolate
#comm = MPI.COMM_WORLD
#mpi_size = comm.Get_size()
#mpi_rank = comm.Get_rank()

status_local = np.zeros(npts_interp)
misloc_local = np.zeros(npts_interp)

status_interp = np.zeros(npts_interp)
misloc_interp = np.zeros(npts_interp)
model_interp = np.zeros((nmodel,npts_interp))

status_local[:] = -1
misloc_local[:] = np.inf

status_interp[:] = -1
misloc_interp[:] = np.inf
model_interp[:] = np.nan

# GLL
xigll, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
yigll, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
zigll, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

#--- loop over each SEM mesh
#for iproc in range(mpi_rank,nproc,mpi_size):
for iproc in range(nproc):

  print("====== proc# ", iproc)
  sys.stdout.flush()

  #--- read in target SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir, iproc)
  mesh_data = sem_mesh_read(mesh_file)

  #--- read model values of the contributing mesh slice
  gll_dims = mesh_data['gll_dims']
  model_gll = np.zeros((nmodel,)+gll_dims)
  for imodel in range(nmodel):
    model_tag = model_names[imodel]
    model_file = "%s/proc%06d_reg1_%s.bin"%(model_dir, iproc, model_tag)
    with FortranFile(model_file, 'r') as f:
      # note: must use fortran convention when reshape to N-D array!!!
      model_gll[imodel,:,:,:,:] = np.reshape(f.read_ints(dtype='f4'), gll_dims, order='F')

  #--- locate each points
  status_local, ispec_local, uvw_local, misloc_local, misratio_local = sem_locate_points_hex27(mesh_data, xyz_interp)

  # merge interpolation results from the current mesh slice into the final results based on misloc and status  
  # index selection for merge: (not located inside an element yet) and (located for the current mesh slice) and (smaller misloc or located inside an element in this mesh slice)
  ii = (status_interp != 1) & (status_local != -1) & ( (misloc_local < misloc_interp) | (status_local == 1) )

  status_interp[ii] = status_local[ii]
  misloc_interp[ii] = misloc_local[ii]
  #misratio_interp[ii] = misratio_local[ii]

  ipoint_select = np.nonzero(ii)[0]
  for ipoint in ipoint_select:
    hlagx = lagrange_poly(xigll, uvw_local[0,ipoint])[:,0]
    hlagy = lagrange_poly(yigll, uvw_local[1,ipoint])[:,0]
    hlagz = lagrange_poly(zigll, uvw_local[2,ipoint])[:,0]
    model_interp[:,ipoint] = np.sum(model_gll[:,:,:,:,ispec_local[ipoint]]*hlagx[:,None,None]*hlagy[None,:,None]*hlagz[None,None,:], axis=(1,2,3))

#--- synchronize all processes
#comm.Barrier()

##--- gather info to the root process 
#MPI_TAG_misloc = 10
#MPI_TAG_model = 11
#
#if mpi_rank != 0:
#  comm.Send(misloc, dest=0, tag=MPI_TAG_misloc)
#  comm.Send(model_interp, dest=0, tag=MPI_TAG_model)
#
#else:
#  misloc_copy = np.empty(npts_interp)
#  model_interp_copy = np.empty((nmodel,npts_interp))
#  for iproc in range(1,mpi_size):
#    comm.Recv(misloc_copy, source=iproc, tag=MPI_TAG_misloc) 
#    comm.Recv(model_interp_copy, source=iproc, tag=MPI_TAG_model) 
#    for ipoint in range(npts_interp):
#      if misloc_copy[ipoint] < misloc[ipoint]:
#        misloc[ipoint] = misloc_copy[ipoint]
#        model_interp[:,ipoint] = model_interp_copy[:,ipoint]

#--- output interpolated model
outdata = xr.Dataset()
for imodel in range(nmodel):
   model = xr.DataArray(np.transpose(np.reshape(model_interp[imodel,:],(nlon,nlat))), coords={'lon':grd_lon1, 'lat':grd_lat1}, dims=('lat', 'lon'))
   model.attrs = {'unit':model_units[imodel]}
   outdata[model_names[imodel]] = model

outdata.attrs = {
    'description':"sem_slice_sphere",
    'depth_km':depth_km,
    }

outdata.to_netcdf(out_file)

#dataset = Dataset(out_file ,'w',format='NETCDF4_CLASSIC')
#dataset.description = "sem_slice_sphere"
##
#dataset.createDimension('longitude',nlon)
#dataset.createDimension('latitude',nlat)
##
#longitudes = dataset.createVariable('longitude', np.float32, ('longitude',))
#latitudes = dataset.createVariable('latitude', np.float32, ('latitude',))
#longitudes.units = 'degree_east'
#latitudes.units = 'degree_north'
#longitudes[:] = grd_lon1
#latitudes[:] = grd_lat1
## model values
#for imodel in range(nmodel):
#  model = dataset.createVariable(model_names[imodel], np.float32, ('longitude','latitude',))
#  model.units = model_units[imodel]
#  #model.long_name = model_names[imodel]
#  model[:,:] = model_interp[imodel,:].reshape((nlon,nlat))
## location misfit of interpolation points 
#misloc_data = dataset.createVariable('misloc', np.float32, ('longitude','latitude',))
#misloc_data.units = 'meter'
#misloc_data[:,:] = misloc_interp.reshape((nlon,nlat))
##
#dataset.close()
