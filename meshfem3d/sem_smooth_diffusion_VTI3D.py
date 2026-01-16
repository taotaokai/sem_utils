#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""smoothing specfem3D GLL model in the lateral and vertical directions using Gaussian window of spatially varying size

ECEF x,y,z coordinates

"""
import sys
import time
import argparse

import numpy as np
import numba
from mpi4py import MPI
import pyproj

from meshfem3d_constants import R_EARTH, R_EARTH_KM
from meshfem3d_utils import (
    sem_mesh_read,
    sem_mesh_get_vol_gll,
    sem_mesh_mpi_read,
    read_gll_file,
    write_gll_file,
    gll2glob,
    laplacian_ani3D,
    assemble_MPI_scalar,
    get_gll_weights,
)

comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

# ====== parameters
parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int)  # number of slices
parser.add_argument("mesh_dir")  # <mesh_dir>/*_solver_data[_mpi].bin
parser.add_argument("model_dir")  # directory of model files to smooth
parser.add_argument("model_name")  # <model_dir>/proc***_<model_name>.bin
parser.add_argument("ref_model_dir")  # directory of reference velocity model name
parser.add_argument("ref_model_name")  # reference velocity model name
parser.add_argument(
    "min_period", type=float, help="smoothing length = min_period * ref_model"
)  # minimum resolved period
parser.add_argument("out_dir")  # num of processes
parser.add_argument(
    "--horizontal_scale",
    type=float,
    default=1.0,
    help="scaling factor for horizontal smoothing length",
)
parser.add_argument(
    "--vertical_scale",
    type=float,
    default=1.0,
    help="scaling factor for vertical smoothing length",
)
# control parameters for CG solver
parser.add_argument("--nstep", type=int, default=100, help="number of time steps")  # nt
parser.add_argument(
    "--max_iter", type=int, default=100, help="maximum number of iteration"
)
parser.add_argument(
    "--max_tolerance",
    type=float,
    default=1e-5,
    help="relative residual to stop iteration",
)

args = parser.parse_args()
if mpi_rank == 0:
    print(args)

nproc = args.nproc
mesh_dir = args.mesh_dir  # <mesh_dir>/proc******_solver_data[_mpi].bin
model_dir = args.model_dir  # <model_dir>/proc******_<model_name>.bin
model_name = args.model_name  # e.g. vpv,vsv,rho,qmu,qkappa
ref_model_dir = args.ref_model_dir  #
ref_model_name = args.ref_model_name  # e.g. vpv,vsv,rho,qmu,qkappa
min_period = args.min_period
nt = args.nstep
out_dir = args.out_dir
max_iter = args.max_iter
max_tolerance = args.max_tolerance
horizontal_scale = args.horizontal_scale
vertical_scale = args.vertical_scale

# time range [0, 1]
dt = 1.0 / nt

# ====== smooth each target mesh slice
if mpi_rank == 0:
    print(args)

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

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish reading mesh/model, {elapsed_time=}")
    sys.stdout.flush()

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

# read model and get volume averaged value on global nodes
model_gll = read_gll_file(model_dir, model_name, mpi_rank, shape=gll_dims)
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

# smoothing length
ref_model_gll = read_gll_file(ref_model_dir, ref_model_name, mpi_rank, shape=gll_dims)
ref_wavelength_gll = (
    abs(ref_model_gll * min_period) / R_EARTH_KM
)  # resolving wavelength
# smooth_length is defined as the full width at half maximum (FWHM) of gaussian function
# smooth_length = FWHM = 2*sqrt(2*log2(2)) * sigma
ref_sigma_gll = ref_wavelength_gll / (2 * np.sqrt(2 * np.log(2)))  # to sigma
# kappa = sigma^2 / 2
# in local ENU coordinates (easting, northing, up)
Kv_gll = 0.5 * (ref_sigma_gll * vertical_scale) ** 2  # vertical
Kh_gll = 0.5 * (ref_sigma_gll * horizontal_scale) ** 2  # horizontal
# rotate to ECEF coordinates
Kxx_glob = np.zeros(nglob)
Kyy_glob = np.zeros(nglob)
Kzz_glob = np.zeros(nglob)
Kxy_glob = np.zeros(nglob)
Kxz_glob = np.zeros(nglob)
Kyz_glob = np.zeros(nglob)
# get lat/lon of global nodes
xyz_glob = mesh["xyz_glob"]
xyz = R_EARTH * xyz_glob
# convert to lat,lon,alt
ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")  # ECEF to GPS
lat_glob, lon_glob, alt_glob = ecef2gps.transform(
    xyz[:, 0], xyz[:, 1], xyz[:, 2], radians=True
)
# volumetric average of GLL points onto global nodes
Kh_dv_glob = gll2glob(Kh_gll * dv_gll, nglob, ibool)
assemble_MPI_scalar(
    Kh_dv_glob,
    mesh_mpi["num_interfaces"],
    mesh_mpi["max_nibool_interfaces"],
    mesh_mpi["nibool_interfaces"],
    mesh_mpi["ibool_interfaces"],
    mesh_mpi["my_neighbors"],
)
Kh_glob = Kh_dv_glob / dv_glob

Kv_dv_glob = gll2glob(Kv_gll * dv_gll, nglob, ibool)
assemble_MPI_scalar(
    Kv_dv_glob,
    mesh_mpi["num_interfaces"],
    mesh_mpi["max_nibool_interfaces"],
    mesh_mpi["nibool_interfaces"],
    mesh_mpi["ibool_interfaces"],
    mesh_mpi["my_neighbors"],
)
Kv_glob = Kv_dv_glob / dv_glob


@numba.jit(nopython=True, nogil=True)
def rotate_K(Kxx_glob, Kyy_glob, Kzz_glob, Kxy_glob, Kxz_glob, Kyz_glob):
    rotmat = np.zeros((3, 3))
    kl = np.zeros((3, 3))
    for idof in range(nglob):
        # kl in ENU coordinates
        kl[:, :] = 0
        kl[0, 0] = Kh_glob[idof]
        kl[1, 1] = Kh_glob[idof]
        kl[2, 2] = Kv_glob[idof]
        # rotation matrix (ENU to ECEF)
        coslat = np.cos(lat_glob[idof])
        sinlat = np.sin(lat_glob[idof])
        coslon = np.cos(lon_glob[idof])
        sinlon = np.sin(lon_glob[idof])
        # rotmat[i,m] = ECEF_Vi .dot. ENU_Vm (Vi: i-th unit vector)
        rotmat[0, :] = [-sinlon, -sinlat * coslon, coslat * coslon]
        rotmat[1, :] = [coslon, -sinlat * sinlon, coslat * sinlon]
        rotmat[2, :] = [0, coslat, sinlat]
        # rotate from ENU to ECEF: Kij = A_im * kl_mn * A_jn
        for m in range(3):
            for n in range(3):
                Kxx_glob[idof] += rotmat[0, m] * kl[m, n] * rotmat[0, n]
                Kyy_glob[idof] += rotmat[1, m] * kl[m, n] * rotmat[1, n]
                Kzz_glob[idof] += rotmat[2, m] * kl[m, n] * rotmat[2, n]
                Kxy_glob[idof] += rotmat[0, m] * kl[m, n] * rotmat[1, n]
                Kxz_glob[idof] += rotmat[0, m] * kl[m, n] * rotmat[2, n]
                Kyz_glob[idof] += rotmat[1, m] * kl[m, n] * rotmat[2, n]


rotate_K(Kxx_glob, Kyy_glob, Kzz_glob, Kxy_glob, Kxz_glob, Kyz_glob)
comm.Barrier()

# GLL points weights
zgll, wgll, dlag_dzgll = get_gll_weights()

# trapezoidal method
# M * (u_{n+1} - u_n) = dt * K * (u_{n+1} + u_n) / 2
# (M - dt/2 * K) * u_{n+1} = (M + dt/2 * K) * u_n
# (1 - dt/2 * M^{-1} * K) * u_{n+1} = (1 + dt/2 * M^{-1} * K) * u_n


def Kx(x_glob):
    kx_glob = laplacian_ani3D(
        x_glob,
        Kxx_glob,
        Kyy_glob,
        Kzz_glob,
        Kxy_glob,
        Kxz_glob,
        Kyz_glob,
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
    iter = 0
    threshold = rsold * max_tolerance
    while rsold > threshold:
        if iter > max_iter:
            print(
                f"{mpi_rank=}: stop iteration at {max_iter=}, {rsold=}, {pAp=}. "
                f"unconvergence might occur, consider increase number of steps! {nt=}"
            )
            break
        # for i in range(100):
        Ap = Ax(p)
        pAp = comm.allreduce(sum(p * Ap), op=MPI.SUM)
        alpha = rsold / pAp
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = comm.allreduce(sum(r**2), op=MPI.SUM)
        p = r + (rsnew / rsold) * p
        rsold = rsnew
        iter += 1
        # if mpi_rank == 0:
        #     i += 1
        #     print(f"{i=}, {rsold=}, {pAp=}")
        #     sys.stdout.flush()
    return x, iter


if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== start iterative smoothing (diffusion), {elapsed_time=}")
    sys.stdout.flush()

u = u_glob
for it in range(nt):
    # ku = Kx(u)
    # mu = Mx(u)
    u, iter = solve_cg(u)
    maxamp_u = comm.reduce(max(abs(u)), op=MPI.MAX, root=0)
    if mpi_rank == 0:
        elapsed_time = time.time() - tic
        # print(f"{it=}, {max(abs(mu))=}, {max(abs(ku))=}, {max(abs(b))=}, {max(abs(u))=}, {elapsed_time=}")
        print(f"{it=}, max(abs(u))={maxamp_u}, {iter=}, {elapsed_time=}")
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
write_gll_file(out_dir, model_name, mpi_rank, u[ibool])

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish writing out smoothed model files, {elapsed_time=}")
    # print(f"{du_glob=}, {du_glob.dtype=}")

comm.Barrier()
