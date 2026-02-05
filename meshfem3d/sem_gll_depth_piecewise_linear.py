#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
        description="Create piecewise linear depth-dependent model"
    )
    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument("mesh_dir", help="directory with proc*_reg1_solver_data.bin")
    parser.add_argument("--out_dir", default="./", help="output directory of GLL file")
    parser.add_argument(
        "--out_tag",
        default="gauss",
        help="model tag for output GLL file as proc*_reg1_[out_tag].bin",
    )
    parser.add_argument(
        "--depths", required=True, nargs="+", type=float, help="depth grids in km"
    )
    parser.add_argument(
        "--values", required=True, nargs="+", type=float, help="values at each grids"
    )

    return parser.parse_args()


def read_mesh_file(mesh_dir, iproc):
    mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
    mesh_data = sem_mesh_read(mesh_file)
    return mesh_data


def interpolate_points(xyz_glob, depths, values):
    r = np.sqrt(np.sum(xyz_glob ** 2, axis=1))
    d = (1.0 - r) * R_EARTH_KM
    model = np.interp(d, depths, values)
    return model


def main():
    args = parse_arguments()
    if rank_world == 0:
        print(args)

    assert len(args.depths) == len(args.values)
    assert np.all(np.diff(args.depths) > 0)

    for iproc in range(rank_world, args.nproc, size_world):
        mesh_data = read_mesh_file(args.mesh_dir, iproc)
        model = interpolate_points(mesh_data["xyz_glob"], args.depths, args.values)
        ibool = mesh_data["ibool"]
        write_gll_file(args.out_dir, args.out_tag, iproc, model[ibool], overwrite=True)


if __name__ == "__main__":
    main()
