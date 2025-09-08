#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""smoothing specfem3D GLL model in the lateral and vertical directions using Gaussian window of spatially varying size

smooth the model using Gaussian kernel by solving the heat equation

Gauss(r,sigma) = exp(- 0.5 * r^2 / sigma^2) / (2 * PI * sigma^2)^(d/2)   (d = 3 for 3-D, 2 for 2-D)

Heat equation:

Du/Dt - k * D^2(u)/Dx^2 = 0

I.C. u(t=0) = m (model before smoothing)
B.C. n .dot. Du/Dx = 0 (adiabatic boundary condition)

Fundamental solution: P(x, t; k) = exp(-r^2/(4kt)) / (4*PI*k*t)^(d/2) for whole space

for k = 0.5 * sigma^2 and t = 1: P(x,1;k) = Gauss(r,sigma)

Since the solution of the intitial value problem is u(t) = P(x, t) * u(x, t=0) (* is spatial convolution)

, it approximates the smoothing by convolution with a Gaussian window except difference near the boundary:

(e.g. initial value is even extended from the boundary)

"""
import sys
import time
import argparse

import numpy as np

from scipy.io import FortranFile

from mpi4py import MPI

from meshfem3d_constants import *

from meshfem3d_utils import sem_mesh_read, sem_mesh_get_vol_gll, sem_mesh_mpi_read
from meshfem3d_utils import gll2glob, laplacian_iso, assemble_MPI_scalar, get_gll_weights


# ====== parameters
parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int)  # number of slices
parser.add_argument("mesh_dir", help="e.g. DATABASES_MPI/")  # <mesh_dir>/*_solver_data[_mpi].bin
parser.add_argument("model_dir", help="e.g. DATABASES_MPI/")  #
parser.add_argument("model_name", help="e.g. vsv")  # <model_dir>/proc***_<model_name>.bin
parser.add_argument("smooth_length", type=float, help="full width at half maximum (FWHM) of Gaussian kernel in km")  #
parser.add_argument("nstep", type=int)  # nt
parser.add_argument("out_dir")  # num of processes

args = parser.parse_args()

nproc = args.nproc
mesh_dir = args.mesh_dir  # <mesh_dir>/proc******_solver_data[_mpi].bin
model_dir = args.model_dir  # <model_dir>/proc******_<model_name>.bin
model_name = args.model_name  # e.g. vpv,vsv,rho,qmu,qkappa
smooth_length = args.smooth_length
nt = args.nstep
out_dir = args.out_dir

# time range [0, 1]
dt = 1.0 / nt

smooth_length /= R_EARTH_KM  # non-dimensionalize as SEM
# smooth_length is defined as the full width at half maximum (FWHM) of Gaussian function
# smooth_length = FWHM = 2*sqrt(2*log2(2)) * sigma
# for wavelength of smooth_length, the smoothing kernel amplitude is about 2.84% of the peak value at zero wave number
sigma = smooth_length / (2 * np.sqrt(2 * np.log(2)))
# diffsivity coefficient kappa = sigma^2 / 2 
kappa = 0.5 * sigma**2

# ====== smooth each target mesh slice
comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

if mpi_size != nproc:
    raise Exception(f"{mpi_size=} must euqal {nproc=}!")

if mpi_rank == 0:
    tic = time.time()

# --- read in SEM mesh slice
mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, mpi_rank)
mesh = sem_mesh_read(mesh_file)

mesh_file = "%s/proc%06d_reg1_solver_data_mpi.bin" % (mesh_dir, mpi_rank)
mesh_mpi = sem_mesh_mpi_read(mesh_file)

nspec = mesh["nspec"]
nglob = mesh["nglob"]
gll_dims = mesh["gll_dims"]
ibool = mesh["ibool"]

model_file = "%s/proc%06d_reg1_%s.bin" % (model_dir, mpi_rank, model_name)
with FortranFile(model_file, "r") as f:
    model_gll = np.reshape(f.read_ints(dtype="f4"), gll_dims)

# print(f'{mpi_rank=}, {mesh_mpi["my_neighbors"]=}')

comm.Barrier()
if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish reading mesh/model, {elapsed_time=}")
    sys.stdout.flush()

dv_gll = sem_mesh_get_vol_gll(mesh)
dv_glob = gll2glob(
    dv_gll,
    mesh["nglob"],
    mesh["ibool"],
)
assemble_MPI_scalar(
    dv_glob,
    mesh_mpi["num_interfaces"],
    mesh_mpi["max_nibool_interfaces"],
    mesh_mpi["nibool_interfaces"],
    mesh_mpi["ibool_interfaces"],
    mesh_mpi["my_neighbors"],
)

u_dv_glob = gll2glob(
    model_gll * dv_gll,
    nglob,
    ibool,
)
assemble_MPI_scalar(
    u_dv_glob,
    mesh_mpi["num_interfaces"],
    mesh_mpi["max_nibool_interfaces"],
    mesh_mpi["nibool_interfaces"],
    mesh_mpi["ibool_interfaces"],
    mesh_mpi["my_neighbors"],
)

u_glob = u_dv_glob / dv_glob

comm.Barrier()

# GLL points weights
zgll, wgll, dlag_dzgll = get_gll_weights()

# trapzoidal method
# M * (u_{n+1} - u_n) = dt * K * (u_{n+1} + u_n) / 2
# (M - dt/2 * K) * u_{n+1} = (M + dt/2 * K) * u_n
# (1 - dt/2 * M^{-1} * K) * u_{n+1} = (1 + dt/2 * M^{-1} * K) * u_n


def Kx(x_glob):
    kx_glob = laplacian_iso(
        x_glob,
        kappa,
        wgll,
        dlag_dzgll,
        mesh["ibool"],
        mesh["dxsi_dx"],
        mesh["dxsi_dy"],
        mesh["dxsi_dz"],
        mesh["deta_dx"],
        mesh["deta_dy"],
        mesh["deta_dz"],
        mesh["dgam_dx"],
        mesh["dgam_dy"],
        mesh["dgam_dz"],
        mesh["jacobian"],
    )

    assemble_MPI_scalar(
        kx_glob,
        mesh_mpi["num_interfaces"],
        mesh_mpi["max_nibool_interfaces"],
        mesh_mpi["nibool_interfaces"],
        mesh_mpi["ibool_interfaces"],
        mesh_mpi["my_neighbors"],
    )

    return kx_glob


def solve_cg(u):
    def Ax(x_glob):
        return x_glob - 0.5 * dt * Kx(x_glob) / dv_glob

    b = u + 0.5 * dt * Kx(u) / dv_glob
    x = u.copy()
    r = b - Ax(x)
    p = r.copy()
    rsold = comm.allreduce(sum(r**2), op=MPI.SUM)
    # x = np.zeros_like(b)
    i = 0
    threshold = rsold * 1e-4
    while rsold > threshold:
        print(f"{rsold=}")
        # for i in range(100):
        Ap = Ax(p)
        pAp = comm.allreduce(sum(p * Ap), op=MPI.SUM)
        alpha = rsold / pAp
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = comm.allreduce(sum(r**2), op=MPI.SUM)
        p = r + (rsnew / rsold) * p
        rsold = rsnew
        # if mpi_rank == 0:
        #     i += 1
        #     print(f"{i=}, {rsold=}, {pAp=}")
        #     sys.stdout.flush()
    return x


if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== start iterative smoothing (diffusion), {elapsed_time=}")
    sys.stdout.flush()

u = u_glob
for it in range(nt):
    # ku = Kx(u)
    # mu = Mx(u)
    u = solve_cg(u)
    if mpi_rank == 0:
        elapsed_time = time.time() - tic
        # print(f"{it=}, {max(abs(mu))=}, {max(abs(ku))=}, {max(abs(b))=}, {max(abs(u))=}, {elapsed_time=}")
        print(f"{it=}, {max(abs(u))=}, {elapsed_time=}")
        sys.stdout.flush()

comm.Barrier()

# # forward Euler
# # M * (u_{n+1} - u_n) = dt * K * u_n
# # u_{n+1} = u_n + dt * M^{-1} * K * u_n
# for it in range(nt):
#     u_glob += dt * Kx(u_glob) / dv_glob
#     if mpi_rank == 0:
#         elapsed_time = time.time() - tic
#         print(f"{it=}, {max(abs(u_glob))=}, {elapsed_time=}")
#         sys.stdout.flush()

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== after smoothing, {elapsed_time=}")
    # print(f"{du_glob=}, {du_glob.dtype=}")

# write out smoothed model files
out_file = "%s/proc%06d_reg1_%s.bin" % (out_dir, mpi_rank, model_name)
with FortranFile(out_file, "w") as f:
    # f.write_record(np.array(u_glob[ibool], dtype="f4"))
    f.write_record(np.array(u[ibool], dtype="f4"))

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish writing out smoothed model files, {elapsed_time=}")
    # print(f"{du_glob=}, {du_glob.dtype=}")

comm.Barrier()