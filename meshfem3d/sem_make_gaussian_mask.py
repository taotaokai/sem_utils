#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import argparse
import numpy as np
from scipy.io import FortranFile
import numba
from mpi4py import MPI

from meshfem3d_constants import R_EARTH_KM
from meshfem3d_utils import sem_mesh_read

# MPI initialization
comm_world = MPI.COMM_WORLD
size_world = comm_world.Get_size()
rank_world = comm_world.Get_rank()

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create Gaussian mask around source/receiver locations, " \
        "prod(1 - exp(-0.5 * (x - x_i)**2 / sigma_i**2), i=0,n_mask)"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument("mesh_dir", help="directory with proc*_reg1_solver_data.bin")
    parser.add_argument(
        "point_list",
        help="list of x,y,z[,sigma_km], where x,y,z are non-dimensionalized by R_EARTH (e.g. 6371 km), sigma_km is optional",
    )
    parser.add_argument("--sigma_km", type=float, default=10.0, help="one sigma of Gaussian mask in km")
    parser.add_argument("--reverse", action="store_true", help="reverse mask to 1.0 - mask")
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_mask.bin"
    )

    return parser.parse_args()


def read_mesh_file(mesh_dir, iproc):
    mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
    mesh_data = sem_mesh_read(mesh_file)

    return mesh_data


def read_list(point_list, sigma_km):
    import pandas as pd

    # df = pd.read_csv(point_list, header=None, delimiter=r"\s+")
    df = pd.read_csv(point_list)
    npts = df.shape[0]
    xyz_mask = np.zeros((npts, 3), dtype=np.float32)
    sigma_mask = np.zeros(npts, dtype=np.float32)

    assert 'x' in df.columns
    assert 'y' in df.columns
    assert 'z' in df.columns
    xyz_mask[:, 0] = df['x'].to_numpy(dtype=np.float32)
    xyz_mask[:, 1] = df['y'].to_numpy(dtype=np.float32)
    xyz_mask[:, 2] = df['z'].to_numpy(dtype=np.float32)

    if 'sigma_km' in df.columns:
        sigma_mask = df['sigma_km'].to_numpy(dtype=np.float32)
    else:
        sigma_mask[:] = sigma_km
    sigma_mask /= R_EARTH_KM

    return xyz_mask, sigma_mask


@numba.jit(nopython=True, nogil=True)
def make_gaussian_mask(xyz_glob, xyz_mask, sigma_mask):
    nglob = xyz_glob.shape[0]
    nmask = xyz_mask.shape[0]
    weight_mask = np.ones(nglob, dtype=np.float32)

    sigma_squared = sigma_mask**2

    for iglob in range(nglob):
        weight = 1
        for imask in range(nmask):
            weight *= 1.0 - np.exp(
                -0.5
                * sum((xyz_glob[iglob, :] - xyz_mask[imask, :])**2)
                / sigma_squared[imask]
            )
        weight_mask[iglob] = weight

    return weight_mask


# def make_gaussian_mask(xyz_glob, xyz_mask, sigma_mask):
#     nglob = xyz_glob.shape[0]
#     nmask = xyz_mask.shape[0]
#     weight_mask = np.ones(nglob, dtype=np.float32)
#
#     # Precompute squared sigma values
#     sigma_squared = sigma_mask * sigma_mask
#
#     # Vectorized computation using broadcasting
#     for imask in range(nmask):
#         # Compute squared Euclidean distances for all glob points at once
#         diff = xyz_glob - xyz_mask[imask, :]
#         squared_distances = np.sum(diff * diff, axis=1)
#         # Apply Gaussian kernel and accumulate weights
#         weight_mask *= np.exp(-0.5 * squared_distances / sigma_squared[imask])
#
#     return weight_mask


def main():
    args = parse_arguments()

    if size_world != args.nproc:
        raise ValueError(f"{size_world=} != {args.nproc=}")

    xyz_mask, sigma_mask = read_list(args.point_list, args.sigma_km)

    for iproc in range(rank_world, args.nproc, size_world):
        mesh_data = read_mesh_file(args.mesh_dir, iproc)
        tic = time.time()
        weight_mask = make_gaussian_mask(mesh_data["xyz_glob"], xyz_mask, sigma_mask)
        toc = time.time() - tic
        print("iproc=%d, time=%f" % (iproc, toc))

        if args.reverse:
            weight_mask = 1.0 - weight_mask

        ibool = mesh_data["ibool"]
        out_file = "%s/proc%06d_reg1_mask.bin" % (args.out_dir, iproc)
        with FortranFile(out_file, "w") as f:
            f.write_record(
                np.array(weight_mask[ibool], dtype="f4")
            )


if __name__ == "__main__":
    main()
