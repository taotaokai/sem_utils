import sys
import pandas as pd
import numpy as np
import pyvista as pv

# import simplekml
import argparse
import xarray as xr

from meshfem3d_utils import xyz2latlon_deg, geodetic_lat2geocentric_lat

parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
parser.add_argument("xsection_list", help="file of xsection params")
parser.add_argument(
    "-d",
    "--sem_db",
    help="SEM DATABASES_MPI",
    default="DATABASES_MPI",
)
parser.add_argument(
    "-m",
    "--model",
    nargs="+",
    default=["vpv"],
)
parser.add_argument("-o", "--out", default="slices.nc", help="output NETCDF file")
parser.add_argument(
    "-r",
    "--r_range",
    nargs=2,
    metavar=("rmin", "rmax"),
    help="limits along radial direction",
    default=[0.8, 1.0],
    type=float,
)
parser.add_argument(
    "-n",
    "--ngrid",
    nargs=2,
    metavar=("theta", "r"),
    help="number of grids along %(metavar)s",
    default=[100, 50],
    type=int,
)
parser.add_argument("--vtk", action="store_true", help="Output VTK file")

args = parser.parse_args()
print(args)

xsection_params = pd.read_csv(args.xsection_list)

ntheta, nr = args.ngrid
rmin, rmax = args.r_range
r_grid = np.linspace(rmin, rmax, nr)

ds_xsection = xr.Dataset()

for i, params in xsection_params.iterrows():

    lat = np.deg2rad(params["lat"])
    lon = np.deg2rad(params["lon"])
    azimuth = np.deg2rad(params["azimuth"])
    min_theta = np.deg2rad(params["min_theta"])
    max_theta = np.deg2rad(params["max_theta"])

    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat)
    phi = lon
    ve = np.array([-np.sin(phi), np.cos(phi), 0])
    vn = np.array(
        [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
    )
    vr = np.array(
        [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    vx = np.cos(azimuth) * vn + np.sin(azimuth) * ve  # vector along xsection

    theta_width = max_theta - min_theta
    dtheta = theta_width / (ntheta - 1)
    theta_grid = min_theta + dtheta * np.arange(ntheta)

    points = np.zeros((ntheta, nr, 3))
    for itheta, theta in enumerate(theta_grid):
        # theta = min_theta + itheta * dtheta
        v = np.tan(theta) * vx + vr
        v = v / sum(v**2) ** 0.5
        points[itheta, :, :] = r_grid[:, None] * v[None, :]

    if args.vtk:
        # connectivity
        ncells = (ntheta - 1) * (nr - 1)
        connectivity = np.zeros((ncells, 4), dtype=int)
        iy, ix = np.unravel_index(np.arange(ncells), (ntheta - 1, nr - 1))
        ii = 0
        for dx, dy in (0, 0), (1, 0), (1, 1), (0, 1):
            ind = np.ravel_multi_index((iy + dy, ix + dx), (ntheta, nr))
            connectivity[:, ii] = (ix + dx) + nr * (iy + dy)
            ii += 1
        mesh = pv.UnstructuredGrid(
            {pv.CellType.QUAD: connectivity}, points.reshape((-1, 3))
        )
        mesh.save(f"vert_xsection_{i:02d}.vtk")

    #
    da = xr.DataArray(
        points,
        dims=("theta", "radius", "dim3"),
        coords={
            "theta": ("theta", theta_grid, {"unit": "degree"}),
            "radius": (
                "radius",
                r_grid,
                {"unit": "R_EARTH"},
            ),
            "dim3": ("dim3", [0,1,2], {"unit":"xyz"}),
        },
        attrs={"latc": params["lat"], "lonc": params["lon"], "azimuth": params["azimuth"]},
        name=f"vertical_xsection_{i:02d}",
    )
    da.to_netcdf(f"vert_xsection_{i:02d}.nc")
    # ds_xsection[f"vertical_xsection_{i:02d}"] = da

# ds_xsection.to_netcdf("vertical_xsections.netcdf")