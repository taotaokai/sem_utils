#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from mpi4py import MPI

from meshfem3d_utils import (
    sem_mesh_read,
    sem_mesh_get_vol_gll,
    read_gll_file,
    write_gll_file,
)

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perturb model")

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "mesh_dir", help="directory of mesh files, proc*_reg1_solver_data.bin"
    )
    parser.add_argument(
        "prev_dmodel_dir",
        help="directory of pervious model search direction GLL files, proc*_reg1_[dmodel_tag].bin",
    )
    parser.add_argument(
        "prev_kernel_dir",
        help="directory of previous kernel GLL files, proc*_reg1_[in_tag].bin",
    )
    parser.add_argument(
        "curr_kernel_dir",
        help="directory of current kernel GLL files, proc*_reg1_[in_tag].bin",
    )
    parser.add_argument(
        "out_dir", help="output dir for scaled GLL files, proc*_reg1_[dmodel_tag].bin"
    )
    parser.add_argument(
        "--kernel_tags",
        nargs="+",
        default=["vp_kernel", "vs_kernel"],
        help="tag as in prev,curr_kernel_dir/proc*_reg1_[kernel_tag].bin",
    )
    parser.add_argument(
        "--dmodel_tags",
        nargs="+",
        default=["vp_dmodel", "vs_dmodel"],
        help="tag as in prev_dmodel,out_dir/proc*_reg1_[dmodel_tag].bin",
    )
    parser.add_argument(
        "--scaled_amplitude",
        type=float,
        help="value to scale the maximum amplitude of the output GLL file",
        default=0.1,
    )

    return parser.parse_args()


def process(
    nproc,
    mesh_dir,
    prev_dmodel_dir,
    prev_kernel_dir,
    curr_kernel_dir,
    out_dir,
    kernel_tags,
    dmodel_tags,
    scaled_amplitude=0.1,
):
    """Process and scale GLL files to scaled amplitude."""

    # get search direction (d_k+1) by CG method for maximization of objective function:
    #
    #   d_k+1 = g_k+1 - beta * d_k  (if minimization, d_k+1 = -1 * g_k+1 + beta * d_k)
    #
    #   , where d_k+1 is curr_dmodel,
    #           d_k   is prev_dmodel,
    #           g_k+1 is curr_kernel,
    #           g_k   is prev_kernel,
    #           y_k = g_k+1 - g_k,
    #           beta = sum(g_k+1 * y_k) / sum(d_k * y_k)

    # get CG parameter beta
    gy_l = 0
    dy_l = 0
    for iproc in range(mpi_rank, nproc, mpi_size):
        mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")
        mesh_data = sem_mesh_read(mesh_file)
        vol_gll = sem_mesh_get_vol_gll(mesh_data)
        vol_gll = vol_gll.reshape(-1)
        for ker_tag, dm_tag in zip(kernel_tags, dmodel_tags):
            prev_dmodel = read_gll_file(prev_dmodel_dir, dm_tag, iproc)
            prev_kernel = read_gll_file(prev_kernel_dir, ker_tag, iproc)
            curr_kernel = read_gll_file(curr_kernel_dir, ker_tag, iproc)
            dkernel = curr_kernel - prev_kernel  # y_k = g_k+1 - g_k
            gy_l += np.sum(vol_gll * curr_kernel * dkernel)  # sum(g_k+1 * y_k)
            dy_l += np.sum(vol_gll * prev_dmodel * dkernel)  # sum(d_k * y_k)

    gy = mpi_comm.allreduce(gy_l, op=MPI.SUM)
    dy = mpi_comm.allreduce(dy_l, op=MPI.SUM)
    if mpi_comm == 0:
        # for maximization problem this value should be negative since Hessian is negative definite
        print(f"INFO: sum(d_k * (g_k+1 - g_k) ) = {dy}") 
    cg_beta = gy / dy  # beta = sum(g_k+1 * y_k) / sum(d_k * y_k)

    # get maximum amplitude of curr_dmodel
    max_amp_l = 0
    for iproc in range(mpi_rank, nproc, mpi_size):
        for ker_tag, dm_tag in zip(kernel_tags, dmodel_tags):
            prev_dmodel = read_gll_file(prev_dmodel_dir, dm_tag, iproc)
            curr_kernel = read_gll_file(curr_kernel_dir, ker_tag, iproc)
            curr_dmodel = curr_kernel - cg_beta * prev_dmodel
            max_amp_l = max(max_amp_l, np.max(np.abs(curr_dmodel)))

    max_amp = mpi_comm.allreduce(max_amp_l, op=MPI.MAX)
    scale_factor = scaled_amplitude / max_amp

    # write out scaled curr_dmodel
    for iproc in range(mpi_rank, nproc, mpi_size):
        for ker_tag, dm_tag in zip(kernel_tags, dmodel_tags):
            prev_dmodel = read_gll_file(prev_dmodel_dir, dm_tag, iproc)
            curr_kernel = read_gll_file(curr_kernel_dir, ker_tag, iproc)
            curr_dmodel = curr_kernel - cg_beta * prev_dmodel
            write_gll_file(out_dir, dm_tag, iproc, scale_factor * curr_dmodel)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    assert len(args.kernel_tags) == len(args.dmodel_tags)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.mesh_dir,
            args.prev_dmodel_dir,
            args.prev_kernel_dir,
            args.curr_kernel_dir,
            args.out_dir,
            args.kernel_tags,
            args.dmodel_tags,
            scaled_amplitude=args.scaled_amplitude,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
