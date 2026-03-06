#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Smooth SPECFEM3D GLL model in lateral and vertical directions using Gaussian windows.

Uses spatially varying Gaussian windows with diffusion-based smoothing in ECEF coordinates.
"""
import argparse
import sys
import time

import numba
import numpy as np
import pyproj
from mpi4py import MPI

from meshfem3d_constants import R_EARTH, R_EARTH_KM
from meshfem3d_utils import (
    assemble_MPI_scalar,
    gll2glob,
    get_gll_weights,
    laplacian_ani3D,
    read_gll_file,
    sem_mesh_get_vol_gll,
    sem_mesh_mpi_read,
    sem_mesh_read,
    write_gll_file,
)

comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Smooth SPECFEM3D GLL model using diffusion"
    )
    parser.add_argument("nproc", type=int, help="number of processes")
    parser.add_argument("mesh_dir", help="directory containing mesh solver data files")
    parser.add_argument("model_dir", help="directory containing model files to smooth")
    parser.add_argument("model_name", help="model name (e.g., vpv, vsv, rho, qmu, qkappa)")
    parser.add_argument("out_dir", help="output directory for smoothed model files")
    parser.add_argument(
        "--horizontal_length",
        type=float,
        default=1.0,
        help="horizontal smoothing length in km (Gaussian FWHM)",
    )
    parser.add_argument(
        "--vertical_length",
        type=float,
        default=1.0,
        help="vertical smoothing length in km (Gaussian FWHM)",
    )
    parser.add_argument(
        "--nstep",
        type=int,
        default=100,
        help="number of time steps",
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=100,
        help="maximum number of CG iterations",
    )
    parser.add_argument(
        "--max_tolerance",
        type=float,
        default=1e-5,
        help="relative residual tolerance to stop CG iteration",
    )
    return parser.parse_args()


args = parse_arguments()
if mpi_rank == 0:
    print(args)

nproc = args.nproc
mesh_dir = args.mesh_dir
model_dir = args.model_dir
model_name = args.model_name
nt = args.nstep
out_dir = args.out_dir
max_iter = args.max_iter
max_tolerance = args.max_tolerance
horizontal_length = args.horizontal_length
vertical_length = args.vertical_length

dt = 1.0 / nt

if mpi_size != nproc:
    raise ValueError(f"{mpi_size=} must equal {nproc=}!")

if mpi_rank == 0:
    tic = time.time()

mesh_file = f"{mesh_dir}/proc{mpi_rank:06d}_reg1_solver_data.bin"
mesh = sem_mesh_read(mesh_file)

mesh_file = f"{mesh_dir}/proc{mpi_rank:06d}_reg1_solver_data_mpi.bin"
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
u_glob = u_dv_glob / dv_glob

sigma_h = horizontal_length / R_EARTH_KM / (2 * np.sqrt(2 * np.log(2)))
sigma_v = vertical_length / R_EARTH_KM / (2 * np.sqrt(2 * np.log(2)))

Kh_glob = np.full(nglob, 0.5 * sigma_h**2)
Kv_glob = np.full(nglob, 0.5 * sigma_v**2)

Kxx_glob = np.zeros(nglob)
Kyy_glob = np.zeros(nglob)
Kzz_glob = np.zeros(nglob)
Kxy_glob = np.zeros(nglob)
Kxz_glob = np.zeros(nglob)
Kyz_glob = np.zeros(nglob)

xyz_glob = mesh["xyz_glob"]
xyz = R_EARTH * xyz_glob
ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")
lat_glob, lon_glob, _ = ecef2gps.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2], radians=True)


@numba.jit(nopython=True, nogil=True)
def rotate_K(Kxx_glob, Kyy_glob, Kzz_glob, Kxy_glob, Kxz_glob, Kyz_glob):
    """Rotate diffusion tensor from ENU to ECEF coordinates."""
    rotmat = np.zeros((3, 3))
    kl = np.zeros((3, 3))
    for idof in range(nglob):
        kl[:, :] = 0
        kl[0, 0] = Kh_glob[idof]
        kl[1, 1] = Kh_glob[idof]
        kl[2, 2] = Kv_glob[idof]

        coslat = np.cos(lat_glob[idof])
        sinlat = np.sin(lat_glob[idof])
        coslon = np.cos(lon_glob[idof])
        sinlon = np.sin(lon_glob[idof])

        rotmat[0, :] = [-sinlon, -sinlat * coslon, coslat * coslon]
        rotmat[1, :] = [coslon, -sinlat * sinlon, coslat * sinlon]
        rotmat[2, :] = [0, coslat, sinlat]

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

zgll, wgll, dlag_dzgll = get_gll_weights()


def Kx(x_glob):
    """Apply anisotropic Laplacian operator."""
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
    """Solve linear system using conjugate gradient method."""

    def Ax(x_glob):
        return x_glob - 0.5 * dt * Kx(x_glob) / dv_glob

    b = u + 0.5 * dt * Kx(u) / dv_glob
    x = u.copy()
    r = b - Ax(x)
    p = r.copy()
    rsold = comm.allreduce(sum(r**2), op=MPI.SUM)

    niter = 0
    threshold = rsold * max_tolerance
    while rsold > threshold:
        if niter > max_iter:
            print(
                f"{mpi_rank=}: stop iteration at {max_iter=}, {rsold=}, {pAp=}. "
                f"Non-convergence might occur, consider increasing number of steps! {nt=}"
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
        niter += 1
    return x, niter


if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== start iterative smoothing (diffusion), {elapsed_time=}")
    sys.stdout.flush()

u = u_glob
for it in range(nt):
    u, niter = solve_cg(u)
    maxamp_u = comm.reduce(max(abs(u)), op=MPI.MAX, root=0)
    if mpi_rank == 0:
        elapsed_time = time.time() - tic
        print(f"{it=}, max(abs(u))={maxamp_u}, {niter=}, {elapsed_time=}")
        sys.stdout.flush()

comm.Barrier()

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== after smoothing, {elapsed_time=}")

write_gll_file(out_dir, model_name, mpi_rank, u[ibool])

if mpi_rank == 0:
    elapsed_time = time.time() - tic
    print(f"====== finish writing out smoothed model files, {elapsed_time=}")

comm.Barrier()
