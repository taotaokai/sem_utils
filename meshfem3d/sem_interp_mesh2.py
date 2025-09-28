#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""create horizontal slice of SEM model at a given depth"""
import sys
# import warnings
import time
import argparse

import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import interpn
from mpi4py import MPI
# import numba

# from meshfem3d_constants import NGLLX, NGLLY, NGLLZ, GAUSSALPHA, GAUSSBETA
from meshfem3d_constants import (
    IFLAG_CRUST,
    IFLAG_80_MOHO,
    IFLAG_220_80,
    IFLAG_670_220,
    IFLAG_DUMMY,
)

from meshfem3d_utils import (
    # get_gll_weights,
    # lagrange_poly,
    # interp1d_linear,
    # interp_model_gll,
    # sem_mesh_read,
    # sem_mesh_locate_points,
    sem_mesh_interp_model,
)

comm_world = MPI.COMM_WORLD
size_world = comm_world.Get_size()
rank_world = comm_world.Get_rank()

# ====== parameters

parser = argparse.ArgumentParser()

parser.add_argument("nproc_source", help="input mesh nproc", type=int)
parser.add_argument("mesh_dir_source", help="input mesh dir")
parser.add_argument("model_dir_source", help="input model dir")

parser.add_argument("nproc_target", help="nproc", type=int)
parser.add_argument("mesh_dir_target", help="mesh dir for interpolate")
parser.add_argument("model_dir_target", help="output model dir")

parser.add_argument(
    "model_tags", nargs="+", help="model tags to interpolate (e.g. vsv vsh)"
)
parser.add_argument(
    "--method",
    default="linear",
    choices=["gll", "linear"],
    help="Choose interpolation method (gll, linear)",
)
parser.add_argument(
    "--nproc_per_slice",
    default=1,
    help="each target mesh slice is run on nsub processes",
    type=int,
)
parser.add_argument(
    "--output_misloc", action="store_true", help="write out misloc info"
)

args = parser.parse_args()
if rank_world == 0:
    print(args)
    sys.stdout.flush()

nproc_source = args.nproc_source
mesh_dir_source = args.mesh_dir_source
model_dir_source = args.model_dir_source

nproc_target = args.nproc_target
mesh_dir_target = args.mesh_dir_target
model_dir_target = args.model_dir_target

model_tags = args.model_tags

nproc_per_slice = args.nproc_per_slice

# merge regions
idoubling_merge = []
# In SETibet case, since I use a velocity gradient across Moho and no mesh boundary at Moho, treat IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80 as the same region
# idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220]
# idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80]

# ====== interpolate
if size_world % nproc_per_slice != 0:
    raise ValueError(
        f"mpi_size({size_world}) must be divisible by nproc_per_slice({nproc_per_slice})"
    )

my_group = rank_world // nproc_per_slice
num_groups = size_world // nproc_per_slice

comm_group = comm_world.Split(my_group, rank_world)
rank_group = comm_group.Get_rank()

if comm_group != MPI.COMM_NULL:

    for iproc_target in range(my_group, nproc_target, num_groups):
        
        if rank_group == 0:
            tic = time.time()

        sem_mesh_interp_model(
            comm_group,
            iproc_target,
            mesh_dir_target,
            model_dir_target,
            nproc_source,
            mesh_dir_source,
            model_dir_source,
            model_tags,
            idoubling_merge=idoubling_merge,
            method=args.method,
            output_misloc=args.output_misloc,
        )

        if rank_group == 0:
            elapsed_time = time.time() - tic
            print(f"{iproc_target=:03d}, {elapsed_time=:8.3f} seconds")
            sys.stdout.flush()