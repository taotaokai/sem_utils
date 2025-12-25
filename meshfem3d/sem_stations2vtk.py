#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert station file to VTK"""
import pandas as pd
import numpy as np
import pyvista as pv
import pyproj

# import simplekml
import argparse
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

from meshfem3d_constants import R_EARTH

parser = argparse.ArgumentParser()

# parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
parser.add_argument("station_file", help="file of stations")
parser.add_argument("--topo_nc", default=None, help="Topography .nc file")
parser.add_argument("--vtk", default="stations.vtk", help="Output VTK file")

args = parser.parse_args()
print(args)

# ====== read in topo file
# bedrock = Dataset(args.etopo_bedrock_nc, "r", format="NETCDF4")
# geoid = Dataset(args.etopo_geoid_nc, "r", format="NETCDF4")
# etopo_heights = np.array(bedrock.variables["z"]) + np.array(
#     geoid.variables["z"]
# )  # z(lat, lon) height from WGS84 ellipsoid
# etopo_lons = np.array(geoid.variables["lon"])
# etopo_lats = np.array(geoid.variables["lat"])
# etopo_interp = RegularGridInterpolator(
#     (etopo_lats, etopo_lons), etopo_heights, bounds_error=False, fill_value=np.nan
# )
# bedrock.close()
# geoid.close()

# ====== read in station file
df = pd.read_csv(args.station_file, sep=r"\s+", header=None, comment="#")
station_lats = df.iloc[:, 2].values
station_lons = df.iloc[:, 3].values
station_deps = df.iloc[:, 5].values

# ====== calculate station altitudes
if args.topo_nc is not None:
    topo_nc = Dataset(args.topo_nc, "r", format="NETCDF4")
    topo_heights = np.array(topo_nc.variables["z"])
    topo_lons = np.array(topo_nc.variables["lon"])
    topo_lats = np.array(topo_nc.variables["lat"])
    topo_interp = RegularGridInterpolator(
        (topo_lats, topo_lons), topo_heights, bounds_error=False, fill_value=np.nan
    )
    station_topo = topo_interp((station_lats, station_lons))
else:
    # ignore topography
    station_topo = 0.0
station_alts = station_topo - station_deps

# ====== convert to ECEF
ecef2gps = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")  # ECEF to GPS
x, y, z = ecef2gps.transform(station_lats, station_lons, station_alts)

# save to vtk
points = np.vstack((x, y, z)).T
points /= R_EARTH # non-dimensionalized by R_EARTH as done in SEM
mesh = pv.PolyData(points)
mesh.save(args.vtk)