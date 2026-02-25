#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse

import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import RegularGridInterpolator

from netCDF4 import Dataset
import pyproj

from meshfem3d_constants import R_EARTH
from meshfem3d_utils import sem_mesh_read

#====== parameters
parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int)
parser.add_argument("model_file", help="EMC model file (.nc)")
parser.add_argument("--mesh_dir", default="DATABASES_MPI", help="mesh dir")
parser.add_argument("--model_tag", default="VSV", help="parameter name in EMC model")
parser.add_argument("--etopo_dir", default="ETOPO1", help="ETOPO dir")
parser.add_argument("--out_dir", default="interp_model", help="output dir")
parser.add_argument("--out_tag", default="vsv", help="output model name")

args = parser.parse_args()
print(args)

# mesh files for which to get interpolated values
nproc_target = args.nproc
mesh_dir_target = args.mesh_dir

model_file = args.model_file
model_name = args.model_tag
etopo_dir = args.etopo_dir
out_dir = args.out_dir
out_tag = args.out_tag

ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")  # ECEF to GPS

#====== read in topo file
bedrock = Dataset(os.path.join(etopo_dir, "ETOPO_2022_v1_60s_N90W180_bed.nc"), "r", format="NETCDF4")
geoid = Dataset(os.path.join(etopo_dir,"ETOPO_2022_v1_60s_N90W180_geoid.nc"), "r", format="NETCDF4")
etopo1_heights = np.array(bedrock.variables['z']) + np.array(geoid.variables['z']) # z(lat, lon) height from WGS84 ellipsoid
etopo1_lons = np.array(geoid.variables['lon'])
etopo1_lats = np.array(geoid.variables['lat'])
etopo1_interp = RegularGridInterpolator((etopo1_lats, etopo1_lons), etopo1_heights)
bedrock.close()
geoid.close()

#====== read in model file
model = Dataset(model_file, "r")
model_depths = np.array(model.variables['depth'])
model_lats = np.array(model.variables['latitude'])
model_lons = np.array(model.variables['longitude'])
model_values = np.array(model.variables[model_name]) / 1000.0 # m/s to km/s
nlat, nlon = model_lats.size, model_lons.size
model_values = model_values.reshape((-1, nlon, nlat))
model_interp = RegularGridInterpolator((model_depths, model_lons, model_lats), model_values)
model.close()

#====== interpolate
for iproc_target in range(nproc_target):
  print("====== iproc_target ", iproc_target)

  # read in target SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_target, iproc_target)
  mesh_data_target = sem_mesh_read(mesh_file)
  nspec_target = mesh_data_target['nspec']
  ibool_target = mesh_data_target['ibool']
  xyz_glob_target = mesh_data_target['xyz_glob']

  # xyz points to interpolate
  xyz_target = R_EARTH * xyz_glob_target[ibool_target.flatten(), :]

  # convert to lat,lon,alt
  lat, lon, alt = ecef2gps.transform(xyz_target[:,0], xyz_target[:,1], xyz_target[:,2])
  print(np.min(alt), np.max(alt))
  print(np.min(lat), np.max(lat))
  print(np.min(lon), np.max(lon))

  topo = etopo1_interp((lat, lon))
  print(np.min(topo), np.max(topo))

  # convert alt(ellipsoidal height) to depth_km below surface
  depth = -1 * (alt - topo) / 1000.0 # to km

  print(np.min(depth), np.max(depth))

  depth[depth < 0.001] = 0.001

  # interpolate model values
  output_model = model_interp((depth, lon, lat))

  # save interpolated model
  gll_file = "%s/proc%06d_reg1_%s.bin"%(out_dir, iproc_target, out_tag)
  with FortranFile(gll_file, 'w') as f:
    f.write_record(np.array(output_model, dtype='f4'))
