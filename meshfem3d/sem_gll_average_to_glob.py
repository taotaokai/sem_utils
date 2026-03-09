#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Average GLL model values to global points by volume-weighted average of co-located GLL points 

"""
import argparse

from mpi4py import MPI

from meshfem3d_utils import (
    sem_mesh_read,
    sem_mesh_get_vol_gll,
    sem_mesh_mpi_read,
    read_gll_file,
    write_gll_file,
    gll2glob,
    assemble_MPI_scalar,
)

comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

# ====== parameters
parser = argparse.ArgumentParser()

parser.add_argument("--nproc", required=True, type=int)  # number of slices
parser.add_argument("--mesh_dir", required=True)  # <mesh_dir>/*_solver_data[_mpi].bin
parser.add_argument("--model_dir", required=True)  # directory of model files to smooth
parser.add_argument("--model_tags", nargs="+", required=True)  # <model_dir>/proc***_<model_name>.bin
parser.add_argument("--out_dir", required=True)  # <out_dir>/proc***_<model_name>.bin

args = parser.parse_args()
if mpi_rank == 0:
    print(args)

if mpi_size != args.nproc:
    raise Exception(f"{mpi_size=} must euqal {args.nproc=}!")

# --- read in SEM mesh slice
mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (args.mesh_dir, mpi_rank)
mesh = sem_mesh_read(mesh_file)

mesh_file = "%s/proc%06d_reg1_solver_data_mpi.bin" % (args.mesh_dir, mpi_rank)
mesh_mpi = sem_mesh_mpi_read(mesh_file)

nglob = mesh["nglob"]
ibool = mesh["ibool"]
gll_dims = mesh["gll_dims"]

for model_name in args.model_tags:

    model_gll = read_gll_file(args.model_dir, model_name, mpi_rank, shape=gll_dims)

    # get volume averaged value on global nodes
    dv_gll = sem_mesh_get_vol_gll(mesh)
    dv_glob = gll2glob(dv_gll, nglob, ibool)
    assemble_MPI_scalar(
        dv_glob,
        mesh_mpi["num_interfaces"],
        mesh_mpi["max_nibool_interfaces"],
        mesh_mpi["nibool_interfaces"],
        mesh_mpi["ibool_interfaces"],
        mesh_mpi["my_neighbors"],
    )

    u_dv_glob = gll2glob(model_gll * dv_gll, nglob, ibool)
    assemble_MPI_scalar(
        u_dv_glob,
        mesh_mpi["num_interfaces"],
        mesh_mpi["max_nibool_interfaces"],
        mesh_mpi["nibool_interfaces"],
        mesh_mpi["ibool_interfaces"],
        mesh_mpi["my_neighbors"],
    )
    u_glob = u_dv_glob / dv_glob  # volumetric averaged value of GLL points on global nodes

    write_gll_file(args.out_dir, model_name, mpi_rank, u_glob[ibool])