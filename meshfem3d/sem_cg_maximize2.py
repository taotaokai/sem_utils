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
        "prev_model_dir",
        help="directory of pervious model GLL files, proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "curr_model_dir",
        help="directory of current model GLL files, proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "prev_kernel_dir",
        help="directory of previous kernel GLL files, proc*_reg1_[kernel_tag].bin",
    )
    parser.add_argument(
        "curr_kernel_dir",
        help="directory of current kernel GLL files, proc*_reg1_[kernel_tag].bin",
    )
    parser.add_argument(
        "curr_precond_kernel_dir",
        help="directory of current preconditioned kernel GLL files, proc*_reg1_[kernel_tag].bin",
    )
    parser.add_argument(
        "out_dir", help="output dir for scaled GLL files, proc*_reg1_[dmodel_tag].bin"
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=["vp", "vs"],
        help="tag as in prev/curr_model_dir/proc*_reg1_[model].bin",
    )
    parser.add_argument(
        "--kernel_tag",
        default="_kernel",
        help="tag as in prev,curr_kernel_dir/proc*_reg1_[model][kernel_tag].bin",
    )
    parser.add_argument(
        "--dmodel_tag",
        default="_dmodel",
        help="tag as in out_dir/proc*_reg1_[model][dmodel_tag].bin",
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
    prev_model_dir,
    curr_model_dir,
    prev_kernel_dir,
    curr_kernel_dir,
    curr_precond_kernel_dir,
    out_dir,
    model_names,
    kernel_tag,
    dmodel_tag,
    scaled_amplitude=0.1,
):
    """Process and scale GLL files to scaled amplitude."""

    # get search direction (d_k+1) by CG method for maximization of objective function:
    #
    #   d_k+1 = Pg_k+1 - beta * d_k  (if minimization, d_k+1 = -1 * g_k+1 + beta * d_k)
    #
    #   , where 
    #
    #   d_k, d_k+1 is the previous/current search direction,
    #   g_k, g_k+1 is the previous/current gradient,
    #   Pg_k+1 is the current preconditioned gradient,
    #   y_k := g_k+1 - g_k is the gradient difference between two iterations,
    #   beta = sum(Pg_k+1 * y_k) / sum(d_k * y_k)

    # get CG parameter beta
    pgy_l = 0
    dy_l = 0
    for iproc in range(mpi_rank, nproc, mpi_size):
        mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")
        mesh_data = sem_mesh_read(mesh_file)
        vol_gll = sem_mesh_get_vol_gll(mesh_data)
        vol_gll = vol_gll.reshape(-1)
        for model_name in model_names:
            # since the search direction might be clipped by the maximum amplitude,
            # we calculate the acutal model difference
            prev_model = read_gll_file(prev_model_dir, model_name, iproc)
            curr_model = read_gll_file(curr_model_dir, model_name, iproc)
            dmodel = curr_model - prev_model # d_k
            # gradient difference
            ker_name = f"{model_name}{kernel_tag}"
            prev_kernel = read_gll_file(prev_kernel_dir, ker_name, iproc)
            curr_kernel = read_gll_file(curr_kernel_dir, ker_name, iproc)
            dkernel = curr_kernel - prev_kernel  # y_k = g_k+1 - g_k
            # preconditioned gradient
            curr_precond_kernel = read_gll_file(curr_precond_kernel_dir, ker_name, iproc) # Pg_k+1
            pgy_l += np.sum(vol_gll * curr_precond_kernel * dkernel)  # sum(Pg_k+1 * y_k)
            dy_l += np.sum(vol_gll * dmodel * dkernel)  # sum(d_k * y_k)

    pgy = mpi_comm.allreduce(pgy_l, op=MPI.SUM)
    dy = mpi_comm.allreduce(dy_l, op=MPI.SUM)
    if mpi_rank == 0:
        # for maximization problem this value should be negative since Hessian is negative definite
        # and Hessian * d_k ~ g_k+1 - g_k, so sum(d_k * (g_k+1 - g_k) ) = d_k * Hessian * d_k < 0
        print(f"INFO: sum(d_k * (g_k+1 - g_k) ) = {dy} should be negative") 
    cg_beta = pgy / dy  # beta = sum(Pg_k+1 * y_k) / sum(d_k * y_k)

    # get maximum amplitude of curr_dmodel
    max_amp_l = 0
    for iproc in range(mpi_rank, nproc, mpi_size):
        for model_name in model_names:
            prev_model = read_gll_file(prev_model_dir, model_name, iproc)
            curr_model = read_gll_file(curr_model_dir, model_name, iproc)
            dmodel = curr_model - prev_model
            ker_name = f"{model_name}{kernel_tag}"
            curr_precond_kernel = read_gll_file(curr_precond_kernel_dir, ker_name, iproc)
            curr_dmodel = curr_precond_kernel - cg_beta * dmodel
            max_amp_l = max(max_amp_l, np.max(np.abs(curr_dmodel)))

    max_amp = mpi_comm.allreduce(max_amp_l, op=MPI.MAX)
    scale_factor = scaled_amplitude / max_amp

    # write out scaled curr_dmodel
    for iproc in range(mpi_rank, nproc, mpi_size):
        for model_name in model_names:
            prev_model = read_gll_file(prev_model_dir, model_name, iproc)
            curr_model = read_gll_file(curr_model_dir, model_name, iproc)
            dmodel = curr_model - prev_model
            ker_name = f"{model_name}{kernel_tag}"
            curr_precond_kernel = read_gll_file(curr_precond_kernel_dir, ker_name, iproc)
            curr_dmodel = curr_precond_kernel - cg_beta * dmodel
            dmodel_name = f"{model_name}{dmodel_tag}"
            write_gll_file(out_dir, dmodel_name, iproc, scale_factor * curr_dmodel)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.mesh_dir,
            args.prev_model_dir,
            args.curr_model_dir,
            args.prev_kernel_dir,
            args.curr_kernel_dir,
            args.curr_precond_kernel_dir,
            args.out_dir,
            args.models,
            args.kernel_tag,
            args.dmodel_tag,
            scaled_amplitude=args.scaled_amplitude,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
