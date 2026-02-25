#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""create horizontal slice of SEM model at a given depth"""

import argparse
import os

import pandas as pd
import numpy as np

import pyproj

from meshfem3d_utils import (
    read_sem_parfile,
    sem_xieta2vec,
    ecef2latlon_zeroalt,
    sem_mesh_interp_points,
    sem_mesh_get_vol_gll,
    sem_mesh_read,
    read_gll_file,
)

from meshfem3d_constants import R_EARTH


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Get 1-D average model")
    parser.add_argument("--sem_parfile", required=True, help="SEM Par_file")
    parser.add_argument(
        "--mesh_dir",
        required=True,
        help="directory of mesh files, e.g. proc******_solver_data.bin",
    )
    parser.add_argument(
        "--model_dir", required=True, help="directory of model GLL files"
    )
    parser.add_argument(
        "--model_names",
        nargs="+",
        required=True,
        help="tag for GLL file as proc*_reg1_[model_name].bin",
    )
    parser.add_argument(
        "--grid_spacing",
        type=float,
        default=1.0,
        help="surface grid spacing in degrees",
    )
    parser.add_argument(
        "--depth_range",
        nargs=2,
        type=float,
        default=[0.0, 1000.0],
        help="miniumm depth in km",
    )
    parser.add_argument(
        "--depth_interval", type=float, default=5.0, help="depth grid spacing in km"
    )
    parser.add_argument("--out_csv", default="1d_profiles.csv", help="output csv file")

    return parser.parse_args()


def calc_1D_average_model_by_weighted_sum_over_gll_points_in_depth_intervals(
    nproc,
    mesh_dir,
    model_dir,
    model_name,
    out_csv,
    min_depth,
    max_depth,
    depth_interval,
):
    # number of depth intervals
    nz = max(1, (max_depth - min_depth) // depth_interval)
    depth_interp = np.linspace(min_depth, max_depth, nz + 1)
    half_depth_win = (
        np.hstack((depth_interp[1] - depth_interp[0], np.diff(depth_interp))) / 2.0
    )
    ndepth = depth_interp.size

    GPS_ELLPS = "WGS84"
    ecef = pyproj.Proj(proj="geocent", ellps=GPS_ELLPS)
    lla = pyproj.Proj(proj="latlong", ellps=GPS_ELLPS)

    weight_model_sum = np.zeros((nproc, ndepth))
    weight_sum = np.zeros((nproc, ndepth))

    for iproc in range(nproc):
        print("#proc = ", iproc)

        mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
        mesh_data = sem_mesh_read(mesh_file)
        ibool = mesh_data["ibool"]
        vol_gll = sem_mesh_get_vol_gll(mesh_data)

        xyz_glob = mesh_data["xyz_glob"]
        lon, lat, height = pyproj.transform(
            ecef,
            lla,
            xyz_glob[0, :] * R_EARTH,
            xyz_glob[1, :] * R_EARTH,
            xyz_glob[2, :] * R_EARTH,
        )
        depth_gll = -1 * height[ibool] / 1000.0

        model_gll = read_gll_file(model_dir, model_name, iproc)

        for idep in range(ndepth):
            dep = depth_interp[idep]
            min_depth = dep - half_depth_win[idep]
            max_depth = dep + half_depth_win[idep]
            idx = (depth_gll >= min_depth) & (depth_gll <= max_depth)
            weight = 1.0 - np.abs(depth_gll[idx] - dep) / half_depth_win[idep]
            weight_model_sum[iproc, idep] = np.sum(
                weight * vol_gll[idx] * model_gll[idx]
            )
            weight_sum[iproc, idep] = np.sum(weight * vol_gll[idx])

        # out_file = "%s/proc%06d_reg1_%s.txt" % (out_dir, iproc, model_name)
        # with open(out_file, "w") as f:
        #     for idep in range(ndepth):
        #         f.write(
        #             "%12.5e  %12.5e\n" % (depth_interp[idep], model_interp[iproc, idep])
        #         )


def sampling_1D_profiles_over_regular_xi_eta_grid(
    sem_parfile,
    mesh_dir,
    model_dir,
    model_names,
    grid_spacing,
    min_depth,
    max_depth,
    depth_interval,
    out_csv,
):

    # read SEM parameters
    sem_params = read_sem_parfile(sem_parfile)
    mesh_nproc_xi = int(sem_params["NPROC_XI"])
    mesh_nproc_eta = int(sem_params["NPROC_ETA"])
    mesh_central_lat = float(sem_params["CENTER_LATITUDE_IN_DEGREES"].lower().replace("d", "e"))
    mesh_central_lon = float(sem_params["CENTER_LONGITUDE_IN_DEGREES"].lower().replace("d", "e"))
    mesh_width_xi = float(sem_params["ANGULAR_WIDTH_XI_IN_DEGREES"].lower().replace("d", "e"))
    mesh_width_eta = float(sem_params["ANGULAR_WIDTH_ETA_IN_DEGREES"].lower().replace("d", "e"))
    mesh_gamma_rot = float(sem_params["GAMMA_ROTATION_AZIMUTH"].lower().replace("d", "e"))

    # define sampling grids
    # number of blocks along xi/eta
    nxi = max(1, int(np.ceil(mesh_width_xi // grid_spacing)))
    neta = max(1, int(np.ceil(mesh_width_eta // grid_spacing)))
    # block edges along xi/eta
    xi = np.linspace(-mesh_width_xi / 2.0, mesh_width_xi / 2.0, nxi + 1)
    eta = np.linspace(-mesh_width_eta / 2.0, mesh_width_eta / 2.0, neta + 1)
    # block centers along xi/eta
    xi = 0.5 * (xi[:-1] + xi[1:])
    eta = 0.5 * (eta[:-1] + eta[1:])
    # 2-D grid of xi/eta
    grd2_xi, grd2_eta = np.meshgrid(xi, eta)
    # 2-D grid of vectors through xi/eta
    grd2_xyz = sem_xieta2vec(
        mesh_central_lat, mesh_central_lon, mesh_gamma_rot, grd2_xi, grd2_eta
    )
    # 2-D grid of lat/lon
    grd2_lat, grd2_lon = ecef2latlon_zeroalt(
        grd2_xyz[..., 0], grd2_xyz[..., 1], grd2_xyz[..., 2]
    )

    # create 3-D grid by adding depth grids
    nz = max(2, int(np.ceil(abs(max_depth - min_depth) / depth_interval) + 1))
    depths = np.linspace(min_depth, max_depth, nz)  # km -> meters

    grd3_lat = np.zeros((nz, neta, nxi))
    grd3_lon = np.zeros((nz, neta, nxi))
    grd3_dep = np.zeros((nz, neta, nxi))
    grd3_xi = np.zeros((nz, neta, nxi))
    grd3_eta = np.zeros((nz, neta, nxi))

    grd3_xi[:] = grd2_xi
    grd3_eta[:] = grd2_eta
    grd3_lat[:] = grd2_lat
    grd3_lon[:] = grd2_lon
    grd3_dep[:] = depths[:, None, None]

    # get ECEF coordinates
    gps2ecef = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")
    grd3_x, grd3_y, grd3_z = gps2ecef.transform(grd3_lat, grd3_lon, - grd3_dep * 1000.0)

    # interpolate models
    npts = grd3_x.size
    interp_points = np.zeros((npts, 3))
    interp_points[:, 0] = grd3_x.flatten()
    interp_points[:, 1] = grd3_y.flatten()
    interp_points[:, 2] = grd3_z.flatten()
    interp_points /= R_EARTH # non-dimensionalize for SEM mesh

    nproc_tot = mesh_nproc_xi * mesh_nproc_eta
    interp_model, final_status, final_misloc, final_misratio = sem_mesh_interp_points(
        nproc_tot, mesh_dir, model_dir, model_names, interp_points
    )
    nmodel = len(model_names)
    interp_model = np.reshape(interp_model, (nz, neta, nxi, nmodel))

    # save results
    data = {}
    data["xi"] = grd3_xi.flatten()
    data["eta"] = grd3_eta.flatten()
    data["lat"] = grd3_lat.flatten()
    data["lon"] = grd3_lon.flatten()
    data["depth_km"] = grd3_dep.flatten()
    for i, tag in enumerate(model_names):
        data[tag] = interp_model[..., i].flatten()
    data["status"] = final_status.flatten()

    df = pd.DataFrame.from_dict(data)
    df.to_csv(out_csv)
    # get mean 1-D profile
    df_mean = df[df["status"] != -1].groupby("depth_km")[model_names].mean()

    root, ext = os.path.splitext(out_csv)
    df_mean.to_csv(root +  "_mean.csv")


if __name__ == "__main__":
    args = parse_arguments()
    print(args)

    sampling_1D_profiles_over_regular_xi_eta_grid(
        args.sem_parfile,
        args.mesh_dir,
        args.model_dir,
        args.model_names,
        args.grid_spacing,
        args.depth_range[0],
        args.depth_range[1],
        args.depth_interval,
        args.out_csv,
    )