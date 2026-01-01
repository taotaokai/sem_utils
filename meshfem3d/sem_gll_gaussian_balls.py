#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import argparse
import numpy as np
from mpi4py import MPI

from meshfem3d_constants import R_EARTH_KM
from meshfem3d_utils import sem_mesh_read, write_gll_file

# MPI initialization
comm_world = MPI.COMM_WORLD
size_world = comm_world.Get_size()
rank_world = comm_world.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create Gaussian balls at given locations with specified sigma and amplitude"
        "amplitude_i * exp(-0.5 * (x - x_i)**2 / sigma_i**2)"
    )
    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument("mesh_dir", help="directory with proc*_reg1_solver_data.bin")
    parser.add_argument(
        "point_list",
        help="list of x,y,z,[sigma_km,amplitude], where x,y,z are non-dimensionalized by R_EARTH (e.g. 6371 km), sigma_km is optional",
    )
    parser.add_argument("--sigma_km", default=10.0, type=int, help="sigma of Gaussian balls in km")
    parser.add_argument("--amplitude", default=1.0, type=float, help="amplitude of Gaussian balls")
    parser.add_argument("--out_dir", default="./", help="output directory of GLL file")
    parser.add_argument(
        "--out_tag",
        default="gauss",
        help="model tag for output GLL file as proc*_reg1_[out_tag].bin",
    )
    return parser.parse_args()


def read_mesh_file(mesh_dir, iproc):
    mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
    mesh_data = sem_mesh_read(mesh_file)
    return mesh_data


def read_list(point_list, sigma_km, amplitude):
    import pandas as pd

    # df = pd.read_csv(point_list, header=None, delimiter=r"\s+")
    df = pd.read_csv(point_list)
    assert "x" in df.columns
    assert "y" in df.columns
    assert "z" in df.columns
    # assert "sigma_km" in df.columns
    # assert "amplitude" in df.columns

    npts = df.shape[0]
    points = np.zeros((npts, 3))
    points[:, 0] = df["x"].to_numpy(dtype=float)
    points[:, 1] = df["y"].to_numpy(dtype=float)
    points[:, 2] = df["z"].to_numpy(dtype=float)

    if "sigma_km" in df.columns:
        sigmas = df["sigma_km"].to_numpy(dtype=float)
    else:   
        sigmas = np.full(npts, sigma_km)
    sigmas = sigmas / R_EARTH_KM

    if "amplitude" in df.columns:
        amplitudes = df["amplitude"].to_numpy(dtype=float)
    else:
        amplitudes = np.full(npts, amplitude)

    return points, sigmas, amplitudes


# @numba.jit(nopython=True, nogil=True)
def make_gaussian_balls(xyz_glob, points, sigmas, amplitudes):
    nglob = xyz_glob.shape[0]
    model = np.zeros(nglob)

    for i, (sigma, amp) in enumerate(zip(sigmas, amplitudes)):
        model[:] += amp * np.exp(
            -0.5 / sigma**2 * np.sum((xyz_glob[:, :] - points[i, :]) ** 2, axis=-1)
        )

    return model


def main():
    args = parse_arguments()
    if rank_world == 0:
        print(args)

    points, sigmas, amplitudes = read_list(args.point_list, args.sigma_km, args.amplitude)

    for iproc in range(rank_world, args.nproc, size_world):
        mesh_data = read_mesh_file(args.mesh_dir, iproc)
        # tic = time.time()
        model = make_gaussian_balls(mesh_data["xyz_glob"], points, sigmas, amplitudes)
        # toc = time.time() - tic
        # print("iproc=%d, time=%f" % (iproc, toc))
        ibool = mesh_data["ibool"]
        write_gll_file(args.out_dir, iproc, args.out_tag, model[ibool])


if __name__ == "__main__":
    main()
