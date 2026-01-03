#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from scipy.io import FortranFile
from mpi4py import MPI

from meshfem3d_utils import read_gll_file, write_gll_file

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Sum up GLL files from different events with masks"
    )

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "dir_list",
        help="list of directories with proc*_reg1_[model_tag|mask_tag].bin",
    )
    parser.add_argument(
        "model_tag", help="tag for GLL file as proc*_reg1_[model_tag].bin"
    )
    parser.add_argument("out_dir", help="output dir for proc*_reg1_[model_tag].bin")
    parser.add_argument(
        "--ncomp",
        type=int,
        help="number of components in GLL file: array shape [nspec,ngllz,nglly,ngllx,ncomp]",
    )
    parser.add_argument(
        "--mask_tag", default=None, help="tag for GLL file as proc*_reg1_[mask_tag].bin"
    )

    return parser.parse_args()


def read_list(filename):
    with open(filename, "r") as f:
        lines = f.read().splitlines()
    return lines


def process_gll_files(gll_folders, iproc, model_tag, mask_tag=None, ncomp=1):
    """Process and sum all GLL files for a given processor."""
    model_gll = None

    for gll_folder in gll_folders:
        mask_gll = 1.0
        if mask_tag is not None:
            mask_gll = read_gll_file(gll_folder, mask_tag, iproc)
            mask_gll = mask_gll[:, None]  # add an extra last dimension for broadcasting

        data_gll = read_gll_file(gll_folder, model_tag, iproc, shape=(-1, ncomp))

        if model_gll is None:
            model_gll = mask_gll * data_gll
        else:
            model_gll += mask_gll * data_gll

    return model_gll


def main():
    """Main function to orchestrate the GLL file summation process."""
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Read list of directories to process
    gll_folders = read_list(args.dir_list)

    # Distribute work across MPI processes
    for iproc in range(mpi_rank, args.nproc, mpi_size):
        try:
            # Process and sum all GLL files for this processor
            model_gll = process_gll_files(
                gll_folders,
                iproc,
                args.model_tag,
                mask_tag=args.mask_tag,
                ncomp=args.ncomp,
            )

            # Write the summed result
            write_gll_file(args.out_dir, args.model_tag, iproc, model_gll)

        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
