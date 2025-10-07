#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from scipy.io import FortranFile
from mpi4py import MPI

from meshfem3d_utils import sem_mesh_read, sem_mesh_get_vol_gll

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perturb model")

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "in_dir", help="directory of input GLL files, proc*_reg1_[in_tag].bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for scaled GLL files, proc*_reg1_[out_tag].bin"
    )
    parser.add_argument(
        "--in_tags",
        nargs="+",
        default=["vp_kernel", "vs_kernel"],
        help="tag as in in_dir/proc*_reg1_[in_tag].bin",
    )
    parser.add_argument(
        "--out_tags",
        nargs="+",
        default=["vp_dmodel", "vs_dmodel"],
        help="tag as in output_dir/proc*_reg1_[out_tag].bin",
    )
    parser.add_argument(
        "--scaled_amplitude",
        type=float,
        help="maximum amplitude in the output scaled GLL file",
        default=0.1,
    )

    return parser.parse_args()


def read_gll_file(model_dir, model_tag, iproc, dtype="f4"):
    """Read a Fortran unformatted file and return the data."""
    filename = os.path.join(model_dir, f"proc{iproc:06d}_reg1_{model_tag}.bin")
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    with FortranFile(filename, "r") as f:
        data = f.read_reals(dtype=dtype)

    return data


def write_gll_file(model_dir, model_tag, iproc, data, dtype="f4"):
    """Write data to a Fortran unformatted file."""
    filename = os.path.join(model_dir, f"proc{iproc:06d}_reg1_{model_tag}.bin")

    with FortranFile(filename, "w") as f:
        f.write_record(np.array(data, dtype=dtype))


def process(nproc, in_dir, in_tags, out_dir, out_tags, scaled_amplitude=0.1):
    """Process and scale GLL files to scaled amplitude."""

    # get maximum amplitude for all GLL files
    max_amp = -1
    for iproc in range(mpi_rank, nproc, mpi_size):
        for tag in in_tags:
            model_gll = read_gll_file(in_dir, tag, iproc)
            max_amp = max(max_amp, np.max(np.abs(model_gll)))
    max_amp = mpi_comm.allreduce(max_amp, op=MPI.MAX)

    scale_factor = scaled_amplitude / max_amp

    # write out scaled GLL files
    for iproc in range(mpi_rank, nproc, mpi_size):
        for in_tag, out_tag in zip(in_tags, out_tags):
            model_gll = read_gll_file(in_dir, in_tag, iproc)
            write_gll_file(out_dir, out_tag, iproc, scale_factor * model_gll)


def main():
    args = parse_arguments()

    assert len(args.in_tags) == len(args.out_tags)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.in_dir,
            args.in_tags,
            args.out_dir,
            args.out_tags,
            scaled_amplitude=args.scaled_amplitude,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
