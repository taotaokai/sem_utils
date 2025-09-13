import os
import sys
import pandas as pd
import numpy as np
import pyvista as pv
import time
# import simplekml
import argparse
import xarray as xr

from scipy.io import FortranFile

from meshfem3d_utils import geodetic_lat2geocentric_lat

from meshfem3d_utils import (
    get_gll_weights,
    lagrange_poly,
    interp1d_linear,
    sem_mesh_read,
    sem_mesh_locate_points,
)

parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
parser.add_argument("xsection_list", help="file of xsection params")
parser.add_argument(
    "--mesh_dir",
    help="SEM DATABASES_MPI",
    default="DATABASES_MPI",
)
parser.add_argument(
    "--model_dir",
    help="SEM GLL model",
    default="DATABASES_MPI",
)
parser.add_argument(
    "--model_names",
    nargs="+",
    default=["vsv"],
)
parser.add_argument("--nc_file", default="slices.nc", help="output NETCDF file")
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
parser.add_argument("-o", "--out_dir", default="./", help="output directory")
parser.add_argument(
    "method", default='linear', choices=['gll', 'linear'], help="Choose interpolation method (gll, linear)"
)

args = parser.parse_args()
print(args)

nproc = args.nproc
mesh_dir = args.mesh_dir
model_dir = args.model_dir
model_names = args.model_names
nmodel = len(model_names)

ntheta, nr = args.ngrid
rmin, rmax = args.r_range
r_grid = np.linspace(rmin, rmax, nr)

# GLL nodes
zgll, wgll, dlag_dzgll = get_gll_weights()

xsection_params = pd.read_csv(args.xsection_list)
for islice, params in xsection_params.iterrows():
    print(f"{islice=:03d}, {params=}")
    sys.stdout.flush()

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
    tangents = np.zeros((ntheta, 3))
    for itheta, theta in enumerate(theta_grid):
        # theta = min_theta + itheta * dtheta
        v = np.sin(theta) * vx + np.cos(theta) * vr
        # v = v / sum(v**2) ** 0.5
        t = np.cos(theta) * vx - np.sin(theta) * vr
        # t = t / sum(t**2) ** 0.5
        points[itheta, :, :] = r_grid[:, None] * v[None, :]
        tangents[itheta, :] = t

    # array of final results
    points = points.reshape((-1, 3))
    npoints = points.shape[0]
    status_target = np.zeros(npoints, dtype="int")
    status_target[:] = -1
    misloc_target = np.zeros(npoints)
    misloc_target[:] = np.inf
    misratio_target = np.zeros(npoints)
    model_interp = np.zeros((nmodel, npoints))
    model_interp[:] = np.nan

    # -- loop over each slice of source SEM mesh
    for iproc in range(nproc):
        tic = time.time()

        # read in source SEM mesh
        mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")
        mesh_data = sem_mesh_read(mesh_file)

        # read in source model
        gll_dims = mesh_data["gll_dims"]
        model_gll = np.zeros((nmodel,) + gll_dims)
        for i, tag in enumerate(model_names):
            model_file = os.path.join(model_dir, f"proc{iproc:06d}_reg1_{tag}.bin")
            with FortranFile(model_file, "r") as f:
                # here gll_dims = [NSPEC, NGLLZ, NGLLY, NGLLX] which is C convention of
                # Fortran array of [NGLLX, NGLLY, NGLLZ, NSPEC]
                model_gll[i, :, :, :, :] = np.reshape(
                    f.read_ints(dtype="f4"), gll_dims
                )

        # locate target points
        status_all, ispec_all, uvw_all, misloc_all, misratio_all = (
            sem_mesh_locate_points(
                mesh_data,
                points,
                kdtree_num_element=2.0,
                max_dist_ratio=2.0,
            )
        )

        # merge interpolation results of mesh slice (iproc_souce) into
        # the final results based on misloc and status

        # index selection for merge:
        # (not located inside an element yet) and
        # (located for the current mesh slice) and
        # ( smaller misloc or located inside an element in this mesh slice )
        ii = (
            (status_target != 1)
            & (status_all != -1)
            & ((misloc_all < misloc_target) | (status_all == 1))
        )

        status_target[ii] = status_all[ii]
        misloc_target[ii] = misloc_all[ii]
        misratio_target[ii] = misratio_all[ii]

        # slower than index slicing but use less memory
        ipoint_select = np.nonzero(ii)[0]
        for ipoint in ipoint_select:
            # interpolation weights
            if args.method == 'linear':
                wx = interp1d_linear(zgll, uvw_all[ipoint, 0])
                wy = interp1d_linear(zgll, uvw_all[ipoint, 1])
                wz = interp1d_linear(zgll, uvw_all[ipoint, 2])
            elif args.method == 'gll':
                wx = lagrange_poly(zgll, uvw_all[ipoint, 0])
                wy = lagrange_poly(zgll, uvw_all[ipoint, 1])
                wz = lagrange_poly(zgll, uvw_all[ipoint, 2])
            else:
                raise ValueError(f"Unknown interpolation method: {args.method}")
            # get interpolated values
            model_interp[:, ipoint] = np.sum(
                model_gll[:, ispec_all[ipoint], :, :, :]
                * wx[None, None, None, :]
                * wy[None, None, :, None]
                * wz[None, :, None, None],
                axis=(1, 2, 3),
            )

        elapsed_time = time.time() - tic
        print(f"{iproc=:03d}, {elapsed_time=:8.3f} seconds")
        sys.stdout.flush()

    model_interp = model_interp.reshape((nmodel, ntheta, nr))

    points = points.reshape((ntheta, nr, 3))

    if args.vtk:
        for i, tag in enumerate(model_names):
            model = model_interp[i, :, :]
            blocks = pv.MultiBlock()
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
            mesh.point_data[tag] = model.flatten()
            mean_val = np.nanmean(model)
            scale = 2 * 10 * dtheta / mean_val
            # plot 1-D radial profiles
            for itheta in range(10, ntheta, 10):
                m0 = model[itheta, :]
                mask = np.isnan(m0)
                if np.all(mask): # skip profile of all nan's
                    continue
                m = model[itheta, :] - mean_val
                p = points[itheta, :, :]
                t = tangents[itheta, :]
                x = p + scale * m[:, None] * t[None, :]
                line = pv.lines_from_points(x)
                # line.point_data[args.model_name] = m0
                blocks.append(line)
                line = pv.lines_from_points(p)
                line.point_data[tag] = m0
                blocks.append(line)
            blocks.append(mesh)
            blocks.save(f"{args.out_dir}/vert_xsection_{islice:02d}_{tag}.vtmb")

#     ds_xsection = xr.Dataset()
#     da = xr.DataArray(
#         points,
#         dims=("theta", "radius"),
#         coords={
#             "theta": ("theta", theta_grid, {"unit": "degree"}),
#             "radius": (
#                 "radius",
#                 r_grid,
#                 {"unit": "R_EARTH"},
#             ),
#         },
#         attrs={"latc": params["lat"], "lonc": params["lon"], "azimuth": params["azimuth"]},
#         name=f"vertical_xsection_{i:02d}",
#     )
#     da.to_netcdf(f"vert_xsection_{i:02d}.nc")
#     # ds_xsection[f"vertical_xsection_{i:02d}"] = da

# ds_xsection.to_netcdf("vertical_xsections.netcdf")
