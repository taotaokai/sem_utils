#!/usr/bin/env python3
"""
Vertical cross-section slicing tool for SEM mesh data.

This script generates vertical cross-sections through a 3D seismic model
represented as a Spectral Element Method (SEM) mesh. It interpolates model
values along specified cross-sections and outputs the results in NetCDF and/or
VTK formats.
"""

import os
import sys
import time
import argparse
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import xarray as xr
import pyvista as pv
from scipy.io import FortranFile

from meshfem3d_utils import (
    geodetic_lat2geocentric_lat,
    get_gll_weights,
    interp_model_gll,
    sem_mesh_read,
    sem_mesh_locate_points,
    read_gll_file,
)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate vertical cross-sections from SEM mesh data"
    )

    parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
    parser.add_argument(
        "xsection_list",
        help="csv file of xsection params (lat,lon,azimuth,min_theta,max_theta)",
    )

    parser.add_argument(
        "--mesh_dir",
        help="SEM DATABASES_MPI directory",
        default="DATABASES_MPI",
    )
    parser.add_argument(
        "--model_dir",
        help="SEM GLL model directory",
        default="DATABASES_MPI",
    )
    parser.add_argument(
        "--model_names",
        nargs="+",
        default=["vsv"],
        help="Model parameter names to extract",
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
        help="number of grids along %(metavar)s directions",
        default=[100, 50],
        type=int,
    )
    parser.add_argument("--vtk", action="store_true", help="Output VTK files")
    parser.add_argument("-o", "--out_dir", default="./", help="output directory")
    parser.add_argument(
        "--method",
        default="linear",
        choices=["gll", "linear"],
        help="Interpolation method (gll or linear)",
    )

    return parser.parse_args()


def create_cross_section_points(
    params: pd.Series, ntheta: int, radial_grid: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create 3D points for a vertical cross-section.

    Args:
        params: Cross-section parameters (lat, lon, azimuth, min_theta, max_theta)
        ntheta: Number of points in theta direction
        radial_grid: Radial coordinate grid

    Returns:
        Tuple of (points, tangents) arrays
    """
    assert radial_grid.ndim == 1

    # Convert to radians
    lat = np.deg2rad(params["lat"])
    lon = np.deg2rad(params["lon"])
    azimuth = np.deg2rad(params["azimuth"])
    min_theta = np.deg2rad(params["min_theta"])
    max_theta = np.deg2rad(params["max_theta"])

    # Convert geographic to geocentric coordinates
    theta_center = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat)
    phi = lon

    # Unit vectors in spherical coordinate system
    ve = np.array([-np.sin(phi), np.cos(phi), 0])  # East
    vn = np.array(
        [
            -np.cos(theta_center) * np.cos(phi),
            -np.cos(theta_center) * np.sin(phi),
            np.sin(theta_center),
        ]
    )  # North
    vr = np.array(
        [
            np.sin(theta_center) * np.cos(phi),
            np.sin(theta_center) * np.sin(phi),
            np.cos(theta_center),
        ]
    )  # Radial

    # Cross-section direction vector
    vx = np.cos(azimuth) * vn + np.sin(azimuth) * ve

    # Create theta grid
    theta_width = max_theta - min_theta
    dtheta = theta_width / (ntheta - 1) if ntheta > 1 else 0
    theta_grid = min_theta + dtheta * np.arange(ntheta)

    # Generate 3D points
    nr = radial_grid.size
    points = np.zeros((ntheta, nr, 3))
    tangents = np.zeros((ntheta, 3))

    for itheta, theta in enumerate(theta_grid):
        # Position vector
        v = np.sin(theta) * vx + np.cos(theta) * vr
        # Tangent vector
        t = np.cos(theta) * vx - np.sin(theta) * vr

        points[itheta, :, :] = radial_grid[:, None] * v[None, :]
        tangents[itheta, :] = t

    return points, tangents, theta_grid


def interpolate_cross_section(
    nproc: int,
    mesh_dir: str,
    model_dir: str,
    model_names: List[str],
    points: np.ndarray,
    zgll: np.ndarray,
    method: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Process a single mesh slice and interpolate model values.

    Returns:
        Tuple of (status, ispec, uvw, misloc, misratio, model_interp)
    """
    npoints = points.shape[0]
    nmodel = len(model_names)

    # Initialize result arrays
    status_final = np.full(npoints, -1, dtype=int)
    misloc_final = np.full(npoints, np.inf)
    misratio_final = np.full(npoints, np.inf)
    model_final = np.full((nmodel, npoints), np.nan)

    # Loop over each slice of source SEM mesh
    for iproc in range(nproc):
        # Read mesh data
        mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")
        mesh_data = sem_mesh_read(mesh_file)

        # Read model data
        gll_dims = mesh_data["gll_dims"]
        nmodel = len(model_names)
        model_gll = np.zeros((nmodel,) + gll_dims)
        for i, tag in enumerate(model_names):
            model_gll[i, :, :, :, :] = read_gll_file(
                model_dir, tag, iproc, region_code="reg1", dtype="f4", gll_dims=gll_dims
            )

        # Locate points in source mesh slice
        status, ispec, uvw, misloc, misratio = sem_mesh_locate_points(
            mesh_data,
            points,
            kdtree_num_element=2.0,
            max_dist_ratio=2.0,
        )

        # merge interpolation results of source mesh slice into
        # the final results based on misloc and status

        # index selection for merge:
        # (not located inside an element yet) and
        # (located for the current mesh slice) and
        # ( smaller misloc or located inside an element in this mesh slice )
        ii = (
            (status_final != 1)
            & (status != -1)
            & ((misloc < misloc_final) | (status == 1))
        )

        status_final[ii] = status[ii]
        misloc_final[ii] = misloc[ii]
        misratio_final[ii] = misratio[ii]

        ipoint_select = np.nonzero(ii)[0]
        interp_model_gll(
            ipoint_select,
            zgll,
            ispec,
            uvw,
            model_gll,
            model_final,
            method=method,
        )

    return model_final


def create_vtk_output(
    islice: int,
    model_names: List[str],
    model_interp: np.ndarray,
    points: np.ndarray,
    tangents: np.ndarray,
    theta_grid: np.ndarray,
    radial_grid: np.ndarray,
    out_dir: str,
    scale_factor: float = 10.0,
) -> None:
    """
    Create VTK output files for visualization.
    """
    ntheta, nr = len(theta_grid), len(radial_grid)
    ncells = (ntheta - 1) * (nr - 1)

    for i, tag in enumerate(model_names):
        model = model_interp[i, :, :]
        blocks = pv.MultiBlock()

        # Create mesh connectivity
        connectivity = np.zeros((ncells, 4), dtype=int)
        iy, ix = np.unravel_index(np.arange(ncells), (ntheta - 1, nr - 1))

        for ii, (dx, dy) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
            connectivity[:, ii] = (ix + dx) + nr * (iy + dy)

        mesh = pv.UnstructuredGrid(
            {pv.CellType.QUAD: connectivity}, points.reshape((-1, 3))
        )
        mesh.point_data[tag] = model.flatten()

        # Calculate scaling for profile visualization
        mean_val = np.nanmean(model)
        if mean_val != 0:
            scale = 2 * scale_factor * (theta_grid[1] - theta_grid[0]) / mean_val
        else:
            scale = scale_factor

        # Plot 1-D radial profiles
        for itheta in range(10, ntheta, 10):
            m0 = model[itheta, :]
            mask = np.isnan(m0)
            if np.all(mask):  # skip profile of all nan's
                continue

            m = model[itheta, :] - mean_val
            p = points[itheta, :, :]
            t = tangents[itheta, :]
            x = p + scale * m[:, None] * t[None, :]

            # Perturbed profile line
            line = pv.lines_from_points(x)
            blocks.append(line)

            # Original profile line
            line = pv.lines_from_points(p)
            line.point_data[tag] = m0
            blocks.append(line)

        blocks.append(mesh)
        vtk_file = f"vert_xsection_{islice:02d}.vtmb"
        blocks.save(os.path.join(out_dir, vtk_file))


def create_netcdf_output(
    xsection_params: pd.DataFrame,
    model_names: List[str],
    all_model_data: Dict[int, np.ndarray],
    all_points: Dict[int, np.ndarray],
    all_theta_grids: Dict[int, np.ndarray],
    radial_grid: np.ndarray,
    filename: str,
) -> None:
    """
    Create NetCDF output with all cross-sections.
    """
    datasets = []

    for islice in all_model_data.keys():
        params = xsection_params.iloc[islice]
        model_data = all_model_data[islice]
        points = all_points[islice]
        theta_grid = all_theta_grids[islice]
        ntheta, nr = len(theta_grid), len(radial_grid)

        # Create dataset for this cross-section
        data_vars = {}
        for i, model_name in enumerate(model_names):
            data_vars[model_name] = (["theta", "radius"], model_data[i, :, :])

        coords = {
            "theta": (["theta"], theta_grid, {"units": "radians"}),
            "radius": (["radius"], radial_grid, {"units": "earth_radii"}),
            "x": (["theta", "radius"], points[:, :, 0]),
            "y": (["theta", "radius"], points[:, :, 1]),
            "z": (["theta", "radius"], points[:, :, 2]),
        }

        attrs = {
            "lat_center": params["lat"],
            "lon_center": params["lon"],
            "azimuth": params["azimuth"],
            "description": f"Vertical cross-section {islice}",
        }

        ds = xr.Dataset(data_vars, coords=coords, attrs=attrs)
        datasets.append(ds)

    # Combine all datasets
    if datasets:
        combined_ds = xr.concat(datasets, dim="section")
        combined_ds.to_netcdf(filename)


def main():
    """Main execution function."""
    args = parse_arguments()

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Grid parameters
    ntheta, nr = args.ngrid
    rmin, rmax = args.r_range
    radial_grid = np.linspace(rmin, rmax, nr)

    # GLL nodes
    zgll, wgll, dlag_dzgll = get_gll_weights()

    # Read cross-section parameters
    xsection_params = pd.read_csv(args.xsection_list)

    # Storage for all results
    all_model_data = {}
    all_points = {}
    all_theta_grids = {}

    # Process each cross-section
    for islice, params in xsection_params.iterrows():
        print(f"Processing cross-section {islice:03d}")
        sys.stdout.flush()

        # Create cross-section points
        points, tangents, theta_grid = create_cross_section_points(
            params, ntheta, radial_grid
        )

        # Reshape for processing
        points_flat = points.reshape((-1, 3))
        model_interp = interpolate_cross_section(
            args.nproc,
            args.mesh_dir,
            args.model_dir,
            args.model_names,
            points_flat,
            zgll,
            args.method,
        )
        # Reshape results
        model_interp = model_interp.reshape((len(args.model_names), ntheta, nr))

        # Store results
        all_model_data[islice] = model_interp
        all_points[islice] = points
        all_theta_grids[islice] = theta_grid

        # Create VTK output if requested
        if args.vtk:
            create_vtk_output(
                islice,
                args.model_names,
                model_interp,
                points,
                tangents,
                theta_grid,
                radial_grid,
                args.out_dir,
            )

    # Create NetCDF output
    create_netcdf_output(
        xsection_params,
        args.model_names,
        all_model_data,
        all_points,
        all_theta_grids,
        radial_grid,
        os.path.join(args.out_dir, args.nc_file),
    )

    print("Processing complete")


if __name__ == "__main__":
    main()
