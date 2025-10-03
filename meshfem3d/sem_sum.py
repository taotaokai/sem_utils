#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from scipy.io import FortranFile
from mpi4py import MPI

# MPI initialization
comm_world = MPI.COMM_WORLD
size_world = comm_world.Get_size()
rank_world = comm_world.Get_rank()

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Sum up GLL files from different events "
    )

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "dir_list",
        help="list of directories with proc*_reg1_[model_tag].bin",
    )
    parser.add_argument("model_tag", help="tag for GLL file as proc*_reg1_[model_tag].bin")
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[model_tag].bin"
    )

    return parser.parse_args()

def read_list(filename):
    with open(filename, 'r') as f:
        lines = f.read().splitlines()
    return lines

def read_fortran_file(filename, dtype="f4"):
    """Read a Fortran unformatted file and return the data."""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    with FortranFile(filename, "r") as f:
        data = f.read_reals(dtype=dtype)
    return data

def write_fortran_file(filename, data, dtype="f4"):
    """Write data to a Fortran unformatted file."""
    with FortranFile(filename, "w") as f:
        f.write_record(np.array(data, dtype=dtype))

def process_gll_files(gll_folders, iproc, model_tag):
    """Process and sum all GLL files for a given processor."""
    model_gll = None
    
    for gll_folder in gll_folders:
        gll_file = os.path.join(gll_folder, f"proc{iproc:06d}_reg1_{model_tag}.bin")
        if model_gll is None:
            model_gll = read_fortran_file(gll_file)
        else:
            model_gll += read_fortran_file(gll_file)
            
    return model_gll

def main():
    """Main function to orchestrate the GLL file summation process."""
    args = parse_arguments()

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Read list of directories to process
    gll_folders = read_list(args.dir_list)

    # Distribute work across MPI processes
    for iproc in range(rank_world, args.nproc, size_world):
        try:
            # Process and sum all GLL files for this processor
            model_gll = process_gll_files(gll_folders, iproc, args.model_tag)
            
            # Write the summed result to output directory
            out_file = os.path.join(args.out_dir, f"proc{iproc:06d}_reg1_{args.model_tag}.bin")
            write_fortran_file(out_file, model_gll)
            
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
