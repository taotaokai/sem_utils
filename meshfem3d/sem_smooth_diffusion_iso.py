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

(e.g. initial value is evenly extended from the boundary in 1-D case)

"""
import sys
import time
import argparse

import numpy as np

from mpi4py import MPI

from meshfem3d_constants import R_EARTH_KM

from meshfem3d_utils import (
    sem_mesh_read,
    sem_mesh_get_vol_gll,
    sem_mesh_mpi_read,
    read_gll_file,
    write_gll_file,
    gll2glob,
    laplacian_iso,
    assemble_MPI_scalar,
    get_gll_weights,
)

comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

# ====== parameters
parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int)  # number of slices
parser.add_argument(
    "mesh_dir", help="e.g. DATABASES_MPI/"
)  # <mesh_dir>/*_solver_data[_mpi].bin
parser.add_argument("model_dir", help="e.g. DATABASES_MPI/")  #
parser.add_argument(
    "model_name", help="e.g. vsv"
)  # <model_dir>/proc***_<model_name>.bin
parser.add_argument(
    "smooth_length",
    type=float,
    help="full width at half maximum (FWHM) of Gaussian kernel in km",
)
parser.add_argument("out_dir")  # num of processes
parser.add_argument(
    "--nstep", type=int, default=100, help="number of time steps between [0, 1]"
)
parser.add_argument(
    "--ratio", type=float, default=1.0, help="ratio of last time step to first time step"
)
# control parameters for CG solver
parser.add_argument(
    "--max_iter", type=int, default=100, help="maximum number of iteration"
)
parser.add_argument(
    "--max_tolerance",
    type=float,
    default=1e-5,
    help="relative residual to stop iteration",
)
parser.add_argument(
    "--method",
    type=str,
    choices=["backward_euler", "crank_nicolson"],
    default="backward_euler",
    help="time integration method",
)
parser.add_argument(
    "--debug",
    action="store_true",
    help="output immediate results for debugging",
)

args = parser.parse_args()
print(args)

nproc = args.nproc
mesh_dir = args.mesh_dir  # <mesh_dir>/proc******_solver_data[_mpi].bin
model_dir = args.model_dir  # <model_dir>/proc******_<model_name>.bin
model_name = args.model_name  # e.g. vpv,vsv,rho,qmu,qkappa
smooth_length = args.smooth_length
nt = args.nstep
out_dir = args.out_dir
max_iter = args.max_iter
max_tolerance = args.max_tolerance

# arithmetically increasing time step length
#
# dt[0] = a
# dt[1] = a + h
# dt[2] = a + 2*h
# ...
# dt[N-1] = a + (N-1)*h
#
# s.t. 
#  sum(dt[i], i=0, N-1) = 1
#  dt[0] = gamma * dt[N-1]
#
# h = 2 / (N * (N - 1)) * (1 - gamma) / (1 + gamma)
# a = 2 * gamma / (1 + gamma) / N
gamma = args.ratio
assert gamma > 0 and gamma <= 1
h = 2 / (nt * (nt - 1)) * (1 - gamma) / (1 + gamma)
a = 2 * gamma / (1 + gamma) / nt
time_steps = a + np.arange(nt) * h

# smoothing length
smooth_length /= R_EARTH_KM  # non-dimensionalize as in SEM
# smooth_length is defined as the full width at half maximum (FWHM) of Gaussian function
# smooth_length = FWHM = 2*sqrt(2*log2(2)) * sigma
# for wavelength of smooth_length, the smoothing kernel amplitude is about 2.84% of the peak value at zero wave number
sigma = smooth_length / (2 * np.sqrt(2 * np.log(2)))
# diffsivity coefficient kappa = sigma^2 / 2
kappa = 0.5 * sigma**2

# ====== smooth each target mesh slice
if mpi_size != nproc:
    raise Exception(f"{mpi_size=} must euqal {nproc=}!")

if mpi_rank == 0:
    tic = time.time()

# --- read in SEM mesh slice
mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, mpi_rank)
mesh = sem_mesh_read(mesh_file)

mesh_file = "%s/proc%06d_reg1_solver_data_mpi.bin" % (mesh_dir, mpi_rank)
mesh_mpi = sem_mesh_mpi_read(mesh_file)

nglob = mesh["nglob"]
ibool = mesh["ibool"]
gll_dims = mesh["gll_dims"]

model_gll = read_gll_file(model_dir, model_name, mpi_rank, shape=gll_dims)

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish reading mesh/model, {elapsed_time=}")
    sys.stdout.flush()

# get volume averaged value on global nodes
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

if args.debug:
    out_file = "%s/proc%06d_reg1_%s.bin" % (out_dir, mpi_rank, model_name)
    write_gll_file(out_dir, "dv", mpi_rank, dv_glob[ibool])
    write_gll_file(out_dir, f"{model_name}_in", mpi_rank, u_glob[ibool])

comm.Barrier()

# GLL points weights
zgll, wgll, dlag_dzgll = get_gll_weights()

#--- time integration method
# 1. Crank-Nicolson
#   M * (u_{n+1} - u_n) = dt * K * (u_{n+1} + u_n) / 2
#   (M - dt/2 * K) * u_{n+1} = (M + dt/2 * K) * u_n
#   (1 - dt/2 * M^{-1} * K) * u_{n+1} = (1 + dt/2 * M^{-1} * K) * u_n
# 2.  backward Euler
#   M * (u_{n+1} - u_n) = dt * K * u_{n+1}
#   (M - dt * K) * u_{n+1} = M * u_n
#   (1 - dt * M^{-1} * K) * u_{n+1} = u_n

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


# solve for each time step
def solve_cg(u, dt):
    def Ax(x_glob):
        if args.method == "backward_euler":
            return x_glob - dt * Kx(x_glob) / dv_glob
        elif args.method == "crank_nicolson":
            return x_glob - 0.5 * dt * Kx(x_glob) / dv_glob
        else:
            raise ValueError(f"{args.method=} is not supported")

    if args.method == "backward_euler":
        b = u
    elif args.method == "crank_nicolson":
        b = u + 0.5 * dt * Kx(u) / dv_glob
    else:
        raise ValueError(f"{args.method=} is not supported")

    x = u.copy()
    r = b - Ax(x)
    p = r.copy()
    rsold = comm.allreduce(sum(r**2), op=MPI.SUM)
    # x = np.zeros_like(b)
    iter = 0
    threshold = rsold * max_tolerance
    pAp = 0.0
    while rsold > threshold:
        if iter > max_iter:
            print(
                f"{mpi_rank=}: stop iteration at {max_iter=}, {rsold=}, {pAp=}. "
                f"unconvergence might occur, consider increase number of steps! {nt=}"
            )
            break
        Ap = Ax(p)
        pAp = comm.allreduce(sum(p * Ap), op=MPI.SUM)
        alpha = rsold / pAp
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = comm.allreduce(sum(r**2), op=MPI.SUM)
        p = r + (rsnew / rsold) * p
        rsold = rsnew
        iter += 1
    return x


if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== start iterative smoothing (diffusion), {elapsed_time=}")
    sys.stdout.flush()

u = u_glob
for it, dt in enumerate(time_steps): 
    max_u = comm.reduce(max(u), op=MPI.MAX, root=0)
    min_u = comm.reduce(min(u), op=MPI.MIN, root=0)
    if mpi_rank == 0:
        elapsed_time = time.time() - tic
        print(f"{it=}, {dt=}, {min_u=}, {max_u=}, {elapsed_time=}")
        sys.stdout.flush()
    u = solve_cg(u, dt)

comm.Barrier()

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== after smoothing, {elapsed_time=}")

# write out smoothed model files
write_gll_file(out_dir, model_name, mpi_rank, u[ibool], overwrite=True)

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish writing out smoothed model files, {elapsed_time=}")

comm.Barrier()