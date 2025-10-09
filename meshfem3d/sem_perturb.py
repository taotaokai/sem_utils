#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from mpi4py import MPI

from meshfem3d_utils import read_gll_file, write_gll_file

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perturb model")

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "model_dir", help="directory of input GLL files, proc*_reg1_[model_tag].bin"
    )
    parser.add_argument(
        "dmodel_dir",
        help="directory of model perturb GLL files, proc*_reg1_[dmodel_tag].bin",
    )
    parser.add_argument(
        "out_dir",
        help="directory of output perturbed model files, proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "--model_tags",
        nargs="+",
        default=["vp", "vs"],
        help="tag as in model_dir/proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "--dmodel_tags",
        nargs="+",
        default=["dvp", "dvs"],
        help="tag as in perturb_dir/proc*_reg1_[dmodel_tag].bin",
    )
    parser.add_argument(
        "--perturb_type",
        choices=["absolute", "relative", "exponential"],
        default="relative",
        help="perturbation type: relative (m * (1 + dm)), absolute (m + dm), exponential (m * exp(dm))",
    )

    return parser.parse_args()


def process(nproc, model_dir, model_tags, dmodel_dir, dmodel_tags, out_dir, perturb_type="relative"):
    """Process and perturb GLL files."""
    for iproc in range(mpi_rank, nproc, mpi_size):
        for model_tag, dmodel_tag in zip(model_tags, dmodel_tags):
            model_gll = read_gll_file(model_dir, model_tag, iproc)
            dm_gll = read_gll_file(dmodel_dir, dmodel_tag, iproc)
            if perturb_type == "relative":
                out_gll = model_gll * (1.0 + dm_gll)
            elif perturb_type == "absolute":
                out_gll = model_gll + dm_gll
            elif perturb_type == "exponential":
                out_gll = model_gll * np.exp(dm_gll)
            else:
                raise ValueError("Invalid perturb_type: %s" % perturb_type)
            write_gll_file(out_dir, model_tag, iproc, out_gll)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    assert len(args.model_tags) == len(args.dmodel_tags)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.model_dir,
            args.model_tags,
            args.dmodel_dir,
            args.dmodel_tags,
            args.out_dir,
            args.perturb_type,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
