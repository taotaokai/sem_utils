#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""create horizontal slice of SEM model at a given depth"""
import argparse

import numpy as np
import pyproj
import xarray as xr
import pyvista as pv

from meshfem3d_constants import R_EARTH

from meshfem3d_utils import (
    geodetic_lat2geocentric_lat,
    ecef2latlon_zeroalt,
    sem_mesh_interp_points,
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate horizontal cross-section from SEM mesh data"
    )
    parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
    parser.add_argument(
        "mesh_dir",
        help="SEM DATABASES_MPI directory",
        default="DATABASES_MPI",
    )
    parser.add_argument(
        "model_dir",
        help="SEM GLL model directory",
        default="DATABASES_MPI",
    )
    parser.add_argument(
        "--model_names",
        nargs="+",
        default=["vsv"],
        help="Model parameter names to extract",
    )
    parser.add_argument(
        "--central_lat", type=float, default=0, help="Mesh central latitude"
    )
    parser.add_argument(
        "--central_lon", type=float, default=0, help="Mesh central longitude"
    )
    parser.add_argument(
        "--rotation_angle", type=float, default=0, help="Mesh rotation angle"
    )
    parser.add_argument(
        "--width_xi", type=float, default=0, help="Mesh angular width in xi"
    )
    parser.add_argument(
        "--width_eta", type=float, default=0, help="Mesh angular width in eta"
    )
    parser.add_argument("--depth", type=float, default=100, help="xsection depth in km")
    parser.add_argument(
        "-n",
        "--ngrid",
        nargs=2,
        metavar=("xi", "eta"),
        help="number of grids along %(metavar)s directions",
        default=[100, 100],
        type=int,
    )
    parser.add_argument("--vtk", default="xsection.vtk", help="output VTK files")
    parser.add_argument("--nc", default="xsection.nc", help="output NETCDF file")
    parser.add_argument(
        "--method",
        default="linear",
        choices=["gll", "linear"],
        help="Interpolation method (gll or linear)",
    )

    return parser.parse_args()


def create_horizontal_xsection_grids(
    xi,  # degrees
    eta,  # degrees
    depth,  # km
    central_lat=0,  # degrees
    central_lon=0,  # degrees
    rotation_angle=0,  # degrees
):
    """
    Create 3D points for a vertical cross-section.

    Args:
        xi: 1-D array of xi coordinates
        eta: 1-D array of eta coordinates
        depth: depth of the horizontal cross-section in km
        central_lat: Central latitude of the mesh chunk
        central_lon: Central longitude of the mesh chunk
        rotation_angle: rotation angle of the mesh chunk

    Returns:
        xyz[n_xi, n_eta, 3]: coordinates of the 2-D grid points
        lats[n_xi, n_eta]: latitudes of the 2-D grid points
        lons[n_xi, n_eta]: longitudes of the 2-D grid points
        alts[n_xi, n_eta]: altitudes of the 2-D grid points
    """
    gps2ecef = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")

    # Convert to radians
    lat0 = np.deg2rad(central_lat)
    lon0 = np.deg2rad(central_lon)
    gamma = np.deg2rad(rotation_angle)

    # Convert geographic to geocentric coordinates
    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat0)
    phi = lon0

    # Unit vectors in spherical coordinate system
    # Easting vector
    ve = np.array([-np.sin(phi), np.cos(phi), 0])  # East
    # Northing vector
    vn = np.array(
        [
            -np.cos(theta) * np.cos(phi),
            -np.cos(theta) * np.sin(phi),
            np.sin(theta),
        ]
    )
    # Radial vector
    vr = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )
    # unit vector along xi at the mesh center
    v_xi = np.cos(gamma) * ve + np.sin(gamma) * vn
    # unit vector along eta at the mesh center
    v_eta = -np.sin(gamma) * ve + np.cos(gamma) * vn

    # 2-D grid
    xi2, eta2 = np.meshgrid(np.deg2rad(xi), np.deg2rad(eta), indexing="ij")

    l_xi = np.tan(xi2)
    l_eta = np.tan(eta2)
    v = (
        l_xi[:, :, None] * v_xi[None, None, :]
        + l_eta[:, :, None] * v_eta[None, None, :]
        + vr[None, None, :]
    )
    v = v / np.sum(v**2, axis=-1, keepdims=True) ** 0.5

    lat2, lon2 = ecef2latlon_zeroalt(v[..., 0], v[..., 1], v[..., 2])
    alt2 = (
        -1 * np.ones_like(lat2) * depth * 1000.0
    )  # depth to negative ellipsoidal height
    lat2 = np.rad2deg(lat2)
    lon2 = np.rad2deg(lon2)
    xx, yy, zz = gps2ecef.transform(lat2, lon2, alt2)
    xyz = np.zeros(xx.shape + (3,))
    xyz[..., 0] = xx
    xyz[..., 1] = yy
    xyz[..., 2] = zz
    xyz = xyz / R_EARTH

    return xyz, lat2, lon2, alt2


def create_vtk_output(
    model_names,
    model_interp,
    xyz,
    out_file,
):
    """
    Create VTK output files for visualization.
    """
    nxi, neta = xyz.shape[:2]
    ncells = (nxi - 1) * (neta - 1)

    # Create mesh connectivity
    connectivity = np.zeros((ncells, 4), dtype=int)
    iy, ix = np.unravel_index(np.arange(ncells), (nxi - 1, neta - 1))

    for ii, (dx, dy) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        connectivity[:, ii] = neta * (iy + dy) + (ix + dx)

    mesh = pv.UnstructuredGrid({pv.CellType.QUAD: connectivity}, xyz.reshape((-1, 3)))

    for i, tag in enumerate(model_names):
        model = model_interp[:, :, i]
        mesh.point_data[tag] = model.flatten()

    mesh.save(out_file)


def create_netcdf_output(
    model_names, model_data, points, xi, eta, lats, lons, alts, out_file, attrs=None
):
    """
    Create NetCDF output with all cross-sections.
    """
    # datasets = []

    # Create dataset for this cross-section
    data_vars = {}
    for i, model_name in enumerate(model_names):
        data_vars[model_name] = (["xi", "eta"], model_data[:, :, i])

    coords = {
        "xi": (["xi"], xi, {"units": "degree"}),
        "eta": (["eta"], eta, {"units": "degree"}),
    }

    data_vars["latitude"] = (["xi", "eta"], lats)
    data_vars["longitude"] = (["xi", "eta"], lons)
    data_vars["altitude"] = (["xi", "eta"], alts)

    data_vars["x"] = (["xi", "eta"], points[:, :, 0])
    data_vars["y"] = (["xi", "eta"], points[:, :, 1])
    data_vars["z"] = (["xi", "eta"], points[:, :, 2])

    ds = xr.Dataset(data_vars, coords=coords, attrs=attrs)
    ds.to_netcdf(out_file)


def main():
    """Main execution function."""
    args = parse_arguments()

    nmodel = len(args.model_names)

    # Grid parameters
    nxi, neta = args.ngrid
    xi_grids = np.linspace(0, args.width_xi, nxi) - 0.5 * args.width_xi
    eta_grids = np.linspace(0, args.width_eta, neta) - 0.5 * args.width_eta

    # Create grids points on the horizontal cross-section
    xyz, lats, lons, alts = create_horizontal_xsection_grids(
        xi_grids,
        eta_grids,
        args.depth,
        central_lat=args.central_lat,
        central_lon=args.central_lon,
        rotation_angle=args.rotation_angle,
    )

    # Reshape for processing
    points_flat = xyz.reshape((-1,3))
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
    model_interp = model_interp.reshape((nxi, neta, nmodel))

    # Create VTK output if requested
    create_vtk_output(
        args.model_names,
        model_interp,
        xyz,
        args.vtk,
    )

    attrs = {
        "central_lat": args.central_lat,
        "central_lon": args.central_lon,
        "rotation_angle": args.rotation_angle,
        "depth": args.depth,
        "width_xi": args.width_xi,
        "width_eta": args.width_eta,
        "description": f"Horizontal cross-section at depth {args.depth} km",
    }
    create_netcdf_output(
        args.model_names,
        model_interp,
        xyz,
        xi_grids,
        eta_grids,
        lats,
        lons,
        alts,
        args.nc,
        attrs=attrs,
    )


if __name__ == "__main__":
    main()