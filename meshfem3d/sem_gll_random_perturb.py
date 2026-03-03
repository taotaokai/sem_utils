#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from mpi4py import MPI

from meshfem3d_utils import sem_mesh_read, read_gll_file, write_gll_file

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()

PERTURB_METHODS = ["absolute", "relative", "exponential"]


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perturb model")

    parser.add_argument(
        "--nproc", required=True, type=int, help="number of model slices"
    )
    parser.add_argument(
        "--mesh_dir",
        required=True,
        help="Directory containing mesh files (proc*_reg1_solver_data.bin)",
    )
    parser.add_argument(
        "--model_dir",
        required=True,
        help="Directory containing input model GLL files (proc*_reg1_[model_tag].bin)",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Directory for output perturbed model files (proc*_reg1_[model_tag].bin)",
    )
    parser.add_argument(
        "--model_tags",
        nargs="+",
        default=["vp", "vs"],
        help="model tags as in model_dir/proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "--amplitudes",
        nargs="+",
        default=[1.0],
        type=float,
        help="amplitudes of random perturbation",
    )
    parser.add_argument(
        "--method",
        choices=PERTURB_METHODS,
        default="absolute",
        help="perturbation methods: absolute (m + dm), relative (m * (1 + dm)), exponential (m * exp(dm))",
    )

    return parser.parse_args()


def apply_perturbation(
    nglob: int,
    ibool: np.ndarray,
    model_gll: np.ndarray,
    amplitude: float,
    method: str,
) -> np.ndarray:
    """Apply perturbation based on specified type."""
    # dm has same value at boundary points shared among elements
    dm_glob = np.random.uniform(-amplitude, amplitude, nglob)
    # dm is NOT continuse across elements
    # dm_gll = np.random.uniform(-amplitude, amplitude, model_gll.shape)
    if method == "absolute":
        model_out = model_gll + dm_glob[ibool]
    elif method == "relative":
        model_out = model_gll * (1.0 + dm_glob[ibool])
    elif method == "exponential":
        model_out = model_gll * np.exp(dm_glob[ibool])
    else:
        raise ValueError(f"Invalid perturbation method: {method}")
    return model_out


def process(
    nproc: int,
    mesh_dir: str,
    model_dir: str,
    out_dir: str,
    model_tags: list,
    amplitudes: list,
    method: str,
) -> None:
    """
    Process and perturb GLL files.

    Args:
        nproc: Number of model slices
        model_dir: Directory with input model files
        out_dir: Output directory
        model_tags: List of model parameter tags
        method: Method to apply perturbation
    """

    for iproc in range(mpi_rank, nproc, mpi_size):
        mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
        mesh = sem_mesh_read(mesh_file)
        nglob = mesh["nglob"]
        ibool = mesh["ibool"]
        gll_dims = mesh["gll_dims"]
        for tag, amp in zip(model_tags, amplitudes):
            model_gll = read_gll_file(model_dir, tag, iproc, shape=gll_dims)
            out_gll = apply_perturbation(nglob, ibool, model_gll, amp, method)
            write_gll_file(out_dir, tag, iproc, out_gll, overwrite=True)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    # Validate arguments
    if not os.path.exists(args.model_dir):
        raise FileNotFoundError(f"model directory not found: {args.model_dir}")

    if len(args.amplitudes) == 1:
        args.amplitudes = args.amplitudes * len(args.model_tags)
    elif len(args.model_tags) != len(args.amplitudes):
        raise ValueError("Number of model tags and amplitudes must be equal")

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.mesh_dir,
            args.model_dir,
            args.out_dir,
            args.model_tags,
            args.amplitudes,
            args.method,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
