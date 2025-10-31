#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert CMTSOLUTION to have ECEF coordinates and tau"""
import os
import argparse
import numpy as np

# from scipy.interpolate import interpn
from scipy.interpolate import RegularGridInterpolator

from netCDF4 import Dataset
from obspy import UTCDateTime

from pyproj import Transformer

# --- user input
parser = argparse.ArgumentParser()

parser.add_argument("cmt_file", help="input CMTSOLUTION file")
parser.add_argument("out_file", help="output CMTSOLUTION file in ECEF")
parser.add_argument(
    "--topo_ncfile",
    default=None,
    help="topo netcdf file, veritcal reference to WGS84 ellipsoid. require lat,lon,z[lat,lon] variables",
)

args = parser.parse_args()
print(args)

cmt_file = args.cmt_file
out_file = args.out_file
topo_ncfile = args.topo_ncfile

transformer = Transformer.from_crs("EPSG:4326", "EPSG:4978")  # WGS84 to ECEF

# --- read in topo file
if topo_ncfile:
    topo_ds = Dataset(topo_ncfile, "r")
    topo_heights = np.array(topo_ds.variables["z"])
    topo_lons = np.array(topo_ds.variables["lon"])
    topo_lats = np.array(topo_ds.variables["lat"])
    topo_interp = RegularGridInterpolator((topo_lats, topo_lons), topo_heights)
    topo_ds.close()

# --- read in Harvard gCMT
with open(cmt_file, "r") as f:
    lines = [x for x in f.readlines() if not (x.startswith("#"))]

header = lines[0].split()
year = header[1]
month = header[2]
day = header[3]
hour = header[4]
minute = header[5]
second = header[6]

lines = [x.split(":") for x in lines]
event_id = lines[1][1].strip()
time_shift = float(lines[2][1])

tau = float(lines[3][1]) / 1.628  # mimic triangle with gaussian
lat = float(lines[4][1])
lon = float(lines[5][1])
dep = float(lines[6][1]) * 1000.0  # meter

# get WGS84 height
topo = 0.0
if topo_ncfile:
    topo = topo_interp((lat, lon))  # get local topography
height = topo - dep  # WGS84 ellipsoidal height of the earthquake
print(f"{cmt_file=},{lon=},{lat=},{topo=},{height=}")

# convert from lla to ECEF(meters)
x, y, z = transformer.transform(lat, lon, height)
# x, y, z = transformer.transform(lat, lon, -dep)

# centroid time: t0
isotime = "{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z".format(
    year, month, day, hour, minute, second
)
t0 = UTCDateTime(isotime) + time_shift
# modify origin time in header line to have centroid time
header[1] = "{:04d}".format(t0.year)
header[2] = "{:02d}".format(t0.month)
header[3] = "{:02d}".format(t0.day)
header[4] = "{:02d}".format(t0.hour)
header[5] = "{:02d}".format(t0.minute)
minisecond = round(t0.microsecond * 1e-3)
header[6] = "{:02d}.{:03d}".format(t0.second, minisecond)
header_line = " ".join(header)

# moment tensor
# 1,2,3 -> r,theta,phi
# harvard cmt use dyn*cm for moment tensor, *1e-7 to N*m
m11 = float(lines[7][1])
m22 = float(lines[8][1])
m33 = float(lines[9][1])
m12 = float(lines[10][1])
m13 = float(lines[11][1])
m23 = float(lines[12][1])
mt_rtp = np.array([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]]) * 1.0e-7

# coordinate transformation matrix (r,theta,phi) to (x,y,z)
r = (x**2 + y**2 + z**2) ** 0.5
theta = np.arccos(z / r)
phi = np.arctan2(y, x)

sthe = np.sin(theta)
cthe = np.cos(theta)
sphi = np.sin(phi)
cphi = np.cos(phi)

a = np.array(
    [
        [sthe * cphi, cthe * cphi, -1.0 * sphi],
        [sthe * sphi, cthe * sphi, cphi],
        [cthe, -1.0 * sthe, 0.0],
    ]
)

mt_xyz = np.dot(np.dot(a, mt_rtp), np.transpose(a))

# write out new CMTSOLUTION_ECEF
with open(out_file, "w") as fp:
    fp.write("%s\n" % (header_line))
    fp.write("%-18s %s\n" % ("event_name:", event_id))
    fp.write("%-18s %+15.8E\n" % ("t0(s):", 0.0))
    fp.write("%-18s %+15.8E\n" % ("tau(s):", tau))
    fp.write("%-18s %+15.8E\n" % ("x(m):", x))
    fp.write("%-18s %+15.8E\n" % ("y(m):", y))
    fp.write("%-18s %+15.8E\n" % ("z(m):", z))
    fp.write("%-18s %+15.8E\n" % ("Mxx(N*m):", mt_xyz[0, 0]))
    fp.write("%-18s %+15.8E\n" % ("Myy(N*m):", mt_xyz[1, 1]))
    fp.write("%-18s %+15.8E\n" % ("Mzz(N*m):", mt_xyz[2, 2]))
    fp.write("%-18s %+15.8E\n" % ("Mxy(N*m):", mt_xyz[0, 1]))
    fp.write("%-18s %+15.8E\n" % ("Mxz(N*m):", mt_xyz[0, 2]))
    fp.write("%-18s %+15.8E\n" % ("Myz(N*m):", mt_xyz[1, 2]))
