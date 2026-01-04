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
import argparse

import numpy as np
import pandas as pd
import xarray as xr
import pyvista as pv
import mpi4py.MPI as MPI

from meshfem3d_utils import (
    R_EARTH_KM,
    geodetic_lat2geocentric_lat,
    sem_mesh_interp_points,
)

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
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
        "--delta",
        nargs=2,
        metavar=("theta", "r"),
        help="grid spacing along %(metavar)s directions as degrees and km",
        default=[1, 50],
        type=float,
    )
    parser.add_argument("--vtk", action="store_true", help="Output VTK files")
    parser.add_argument("-o", "--out_dir", default="./", help="output directory")
    parser.add_argument(
        "--method",
        default="linear",
        choices=["gll", "linear"],
        help="Interpolation method (gll or linear)",
    )
    parser.add_argument(
        "--profile_interval",
        default=10,
        type=int,
        help="interval of vertical profiles along theta grids",
    )
    parser.add_argument(
        "--profile_ref_val",
        default=None,
        nargs="+",
        type=float,
        help="vtk profile reference value",
    )
    parser.add_argument(
        "--profile_norm_amp",
        default=None,
        nargs="+",
        type=float,
        help="vtk profile normalization amplitude",
    )

    return parser.parse_args()


def create_vertical_xsection_grids(
    angles,  # degrees
    radius,
    central_lat=0,  # degrees
    central_lon=0,  # degrees
    azimuth=0,  # degrees
):
    """
    Create 3D points for a vertical cross-section.

    Args:
        angles: 1-D array of angular coordinates
        radius: 1-D array of radial coordinates
        central_lat: Central latitude of the cross-section
        central_lon: Central longitude of the cross-section
        azimuth: Azimuth angle of the cross-section

    Returns:
        points[na, nr, 3]: grid points as 3-D array (angles.size, radius.size, 3)
        tangents[na, 3]: tangent vectors
    """

    # Convert to radians
    lat0 = np.deg2rad(central_lat)
    lon0 = np.deg2rad(central_lon)
    az = np.deg2rad(azimuth)
    angles = np.deg2rad(angles)

    # Convert geographic to geocentric coordinates
    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat0)
    phi = lon0

    # Unit vectors in spherical coordinate system
    ve = np.array([-np.sin(phi), np.cos(phi), 0])  # East
    vn = np.array(
        [
            -np.cos(theta) * np.cos(phi),
            -np.cos(theta) * np.sin(phi),
            np.sin(theta),
        ]
    )  # North
    vr = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )  # Radial

    # tangent vector along the cross-section
    vx = np.cos(az) * vn + np.sin(az) * ve
    # vz = np.sin(az) * vn - np.cos(az) * ve
    # print(f"normal: {vz}")

    # create 2-D mesh grids
    # na = angles.size
    # nr = radius.size
    # points = np.zeros((na, nr, 3))
    # tangents = np.zeros((na, 3))

    vv = np.sin(angles[:, None]) * vx[None, :] + np.cos(angles[:, None]) * vr[None, :]
    points = radius[None, :, None] * vv[:, None, :]

    tangents = (
        np.cos(angles[:, None]) * vx[None, :] - np.sin(angles[:, None]) * vr[None, :]
    )

    # for i, angle in enumerate(angles):
    #     angle = np.deg2rad(angle)
    #     # radial vector
    #     v = np.sin(angle) * vx + np.cos(angle) * vr
    #     # Tangent vector
    #     t = np.cos(angle) * vx - np.sin(angle) * vr
    #     #
    #     points[i, :, :] = radius[:, None] * v[None, :]
    #     tangents[i, :] = t

    return points, tangents


def create_vtk_output(
    model_names,
    model_interp,
    points,
    tangents,
    angles,
    radius,
    out_dir,
    label="",
    scale_factor=1.0,
    angle_interval=10,
    ref_values=None,
    norm_amplitudes=None,
):
    """
    Create VTK output files for visualization.
    """
    na, nr = len(angles), len(radius)
    ncells = (na - 1) * (nr - 1)

    # Create vtk mesh for cross-section
    connectivity = np.zeros((ncells, 4), dtype=int)
    iy, ix = np.unravel_index(np.arange(ncells), (na - 1, nr - 1))
    for ii, (dx, dy) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        connectivity[:, ii] = (ix + dx) + nr * (iy + dy)
    mesh_xsection = pv.UnstructuredGrid(
        {pv.CellType.QUAD: connectivity}, points.reshape((-1, 3)).astype(np.float32)
    )

    # set depth value on mesh grid
    depth_km = np.zeros((na, nr))
    depth_km[:, :] = (1.0 - radius[None, :]) * R_EARTH_KM
    mesh_xsection.point_data["depth"] = depth_km.flatten().astype(np.float32)

    max_radius = np.max(radius)
    dangle = angles[angle_interval] - angles[0]
    profile_interval = np.deg2rad(dangle) * max_radius

    # fn_root, _ = os.path.splitext(out_file)

    # blocks = pv.MultiBlock()
    for i, tag in enumerate(model_names):
        model = model_interp[:, :, i]

        mesh_xsection.point_data[tag] = model.flatten().astype(np.float32)
        profiles = []

        # Calculate scaling for profile visualization
        if ref_values is None:
            profile_ref_val = np.nanmean(model)
        else:
            profile_ref_val = ref_values[i]
        if norm_amplitudes is None:
            profile_norm_amp = np.nanmax(np.abs(model - profile_ref_val))
        else:
            profile_norm_amp = norm_amplitudes[i]
        print(f"{tag}: ref_val={profile_ref_val}, norm_amp={profile_norm_amp}")

        if profile_norm_amp != 0:
            scale = scale_factor * profile_interval / profile_norm_amp
        else:
            scale = scale_factor

        # Plot 1-D radial profiles
        for itheta in range(0, na, angle_interval):
            m0 = model[itheta, :]
            mask = np.isnan(m0)
            if np.all(mask):  # skip profile of all nan's
                continue

            m = model[itheta, :] - profile_ref_val
            p = points[itheta, :, :]
            t = tangents[itheta, :]
            x = p + scale * m[:, None] * t[None, :]

            # Perturbed profile line
            line1 = pv.lines_from_points(x)
            # Original profile line
            line2 = pv.lines_from_points(p)
            # line2.point_data[tag] = m0

            profiles.extend([line1, line2])
        mesh_profile = pv.merge(profiles)
        out_file = os.path.join(out_dir, f"{tag}_profile_{label}.vtk")
        mesh_profile.save(out_file)

    out_file = os.path.join(out_dir, f"vslice_{label}.vtk")
    mesh_xsection.save(out_file)


def create_netcdf_output(
    model_names, model_data, points, angles, radius, out_file, attrs=None
):
    """
    Create NetCDF output with all cross-sections.
    """
    # datasets = []

    # Create dataset for this cross-section
    data_vars = {}
    for i, model_name in enumerate(model_names):
        data_vars[model_name] = (
            ["theta", "radius"],
            model_data[:, :, i].astype(np.float32),
        )

    coords = {
        "theta": (["theta"], angles, {"units": "degree"}),
        "radius": (["radius"], radius, {"units": f"{R_EARTH_KM} km"}),
    }

    data_vars["x"] = (["theta", "radius"], points[:, :, 0].astype(np.float32))
    data_vars["y"] = (["theta", "radius"], points[:, :, 1].astype(np.float32))
    data_vars["z"] = (["theta", "radius"], points[:, :, 2].astype(np.float32))

    ds = xr.Dataset(data_vars, coords=coords, attrs=attrs)
    ds.to_netcdf(out_file)

    # datasets.append(ds)
    # # Combine all datasets
    # if datasets:
    #     combined_ds = xr.concat(datasets, dim="section")
    #     combined_ds.to_netcdf(filename)


def main():
    """Main execution function."""
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    # Create output directory if it doesn't exist
    if mpi_rank == 0:
        os.makedirs(args.out_dir, exist_ok=True)
    mpi_comm.barrier()

    nmodel = len(args.model_names)

    profile_ref_val = args.profile_ref_val
    profile_norm_amp = args.profile_norm_amp

    if profile_ref_val is not None:
        if len(profile_ref_val) == 1 and nmodel > 1:
            profile_ref_val = profile_ref_val * nmodel
        assert len(profile_ref_val) == nmodel
    if profile_norm_amp is not None:
        if len(profile_norm_amp) == 1 and nmodel > 1:
            profile_norm_amp = profile_norm_amp * nmodel
        assert len(profile_norm_amp) == nmodel

    # Grid parameters
    dtheta, dr = args.delta
    rmin, rmax = args.r_range
    dr = dr / R_EARTH_KM
    nr = int(np.ceil((rmax - rmin) / dr)) + 1
    radius = np.linspace(rmin, rmax, nr)

    # Read cross-section parameters
    xsection_params = pd.read_csv(args.xsection_list, comment="#")

    # Process each cross-section

    # for islice, params in xsection_params.iterrows():
    for islice in range(mpi_rank, len(xsection_params), mpi_size):
        print(f"Processing cross-section {islice:03d}")

        params = xsection_params.iloc[islice]

        lat0 = params["lat"]
        lon0 = params["lon"]
        azimuth = params["azimuth"]

        min_theta = params["min_theta"]
        max_theta = params["max_theta"]
        na = int(np.ceil(max_theta - min_theta) / dtheta) + 1
        print(f"{nr=}, {na=}")
        sys.stdout.flush()
        angles = np.linspace(min_theta, max_theta, na)

        # Create cross-section points
        points, tangents = create_vertical_xsection_grids(
            angles, radius, central_lat=lat0, central_lon=lon0, azimuth=azimuth
        )

        # Reshape for processing
        points_flat = points.reshape((-1, 3))
        model_interp, *_ = sem_mesh_interp_points(
            args.nproc,
            args.mesh_dir,
            args.model_dir,
            args.model_names,
            points_flat,
            idoubling=None,
            method=args.method,
        )
        # Reshape results
        model_interp = model_interp.reshape((na, nr, nmodel))

        # # Store results
        # all_model_data[islice] = model_interp
        # all_points[islice] = points
        # all_theta_grids[islice] = theta_grid

        # Create VTK output if requested
        if args.vtk:
            # out_file = os.path.join(args.out_dir, f"xsection_{islice:04d}.vtk")
            create_vtk_output(
                args.model_names,
                model_interp,
                points,
                tangents,
                angles,
                radius,
                args.out_dir,
                label=f"{islice:04d}",
                angle_interval=args.profile_interval,
                ref_values=profile_ref_val,
                norm_amplitudes=profile_norm_amp,
            )

        # Create NetCDF output
        out_file = os.path.join(args.out_dir, f"xsection_{islice:04d}.nc")
        attrs = {
            "central_lat": lat0,
            "central_lon": lon0,
            "azimuth": azimuth,
            "description": "Vertical cross-section",
        }
        create_netcdf_output(
            args.model_names,
            model_interp,
            points,
            angles,
            radius,
            out_file,
            attrs=attrs,
        )


if __name__ == "__main__":
    main()
