#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import time
import argparse

import numpy as np
from mpi4py import MPI

# from meshfem3d_constants import (
#     IFLAG_CRUST,
#     IFLAG_80_MOHO,
#     IFLAG_220_80,
#     IFLAG_670_220,
#     IFLAG_DUMMY,
# )

from meshfem3d_utils import sem_mesh_interp_mesh

# MPI initialization
comm_world = MPI.COMM_WORLD
size_world = comm_world.Get_size()
rank_world = comm_world.Get_rank()

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create horizontal slice of SEM model at a given depth")
    
    # Input parameters
    parser.add_argument("nproc_source", help="input mesh nproc", type=int)
    parser.add_argument("mesh_dir_source", help="input mesh dir")
    parser.add_argument("model_dir_source", help="input model dir")
    
    # Output parameters
    parser.add_argument("nproc_target", help="nproc", type=int)
    parser.add_argument("mesh_dir_target", help="mesh dir for interpolate")
    parser.add_argument("model_dir_target", help="output model dir")
    
    # Model parameters
    parser.add_argument("model_tags", nargs="+", help="model tags to interpolate (e.g. vsv vsh)")
    
    # Options
    parser.add_argument("--method", default="linear", choices=["gll", "linear"], 
                       help="Choose interpolation method (gll, linear)")
    parser.add_argument("--nproc_per_slice", default=1, type=int,
                       help="each target mesh slice is run on nsub processes")
    parser.add_argument("--output_misloc", action="store_true", 
                       help="write out misloc info")
    
    return parser.parse_args()

def create_mpi_groups(nproc_per_slice):
    if size_world % nproc_per_slice != 0:
        raise ValueError(
            f"MPI size ({size_world}) must be divisible by nproc_per_slice ({nproc_per_slice})"
        )
    
    my_group = rank_world // nproc_per_slice
    num_groups = size_world // nproc_per_slice
    
    comm_group = comm_world.Split(my_group, rank_world)
    return comm_group, my_group, num_groups

def process_mesh_slices(comm_group, my_group, num_groups, args):
    rank_group = comm_group.Get_rank()
    
    # Parameters from arguments
    nproc_target = args.nproc_target
    mesh_dir_target = args.mesh_dir_target
    model_dir_target = args.model_dir_target
    nproc_source = args.nproc_source
    mesh_dir_source = args.mesh_dir_source
    model_dir_source = args.model_dir_source
    model_tags = args.model_tags
    
    # Merge regions configuration
    idoubling_merge = []  # Can be customized as needed
    
    for iproc_target in range(my_group, nproc_target, num_groups):
        if rank_group == 0:
            start_time = time.time()
        
        sem_mesh_interp_mesh(
            comm_group,
            iproc_target,
            mesh_dir_target,
            model_dir_target,
            nproc_source,
            mesh_dir_source,
            model_dir_source,
            model_tags,
            # idoubling_merge=idoubling_merge,
            method=args.method,
            output_misloc=args.output_misloc,
        )
        
        if rank_group == 0:
            elapsed_time = time.time() - start_time
            print(f"Processed slice {iproc_target:03d} in {elapsed_time:8.3f} seconds")

def main():
    args = parse_arguments()
    
    if rank_world == 0:
        print(args)
        sys.stdout.flush()
    
    comm_group, my_group, num_groups = create_mpi_groups(args.nproc_per_slice)
    
    if comm_group != MPI.COMM_NULL:
        process_mesh_slices(comm_group, my_group, num_groups, args)

if __name__ == "__main__":
    main()