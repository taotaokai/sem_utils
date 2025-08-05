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
import numba_mpi

from scipy.io import FortranFile

from mpi4py import MPI

from meshfem3d_utils import sem_mesh_read, sem_mesh_get_vol_gll, sem_mesh_mpi_read
from meshfem3d_constants import *

import gll_library


# @numba.jit
def assemble_MPI_scalar(
    array_glob,
    num_interfaces,
    max_nibool_interfaces,
    nibool_interfaces,
    ibool_interfaces,
    my_neighbors,
):

    # integer, intent(in) :: num_interfaces,max_nibool_interfaces
    # integer, dimension(num_interfaces), intent(in) :: nibool_interfaces,my_neighbors
    # integer, dimension(max_nibool_interfaces,num_interfaces), intent(in) :: ibool_interfaces

    # comm = MPI.COMM_WORLD
    # mpi_rank = comm.Get_rank()
    mpi_rank = numba_mpi.rank()

    buffer_send_scalar = np.zeros(
        (num_interfaces, max_nibool_interfaces), dtype=array_glob.dtype
    )
    buffer_recv_scalar = np.zeros(
        (num_interfaces, max_nibool_interfaces), dtype=array_glob.dtype
    )

    for iinterface in range(num_interfaces):
        npts = nibool_interfaces[iinterface]
        buffer_send_scalar[iinterface, :npts] = array_glob[
            ibool_interfaces[iinterface, :npts]
        ]

    # send messages
    # send_requests = np.zeros((num_interfaces,), dtype=numba_mpi.RequestType)
    # recv_requests = np.zeros((num_interfaces,), dtype=numba_mpi.RequestType)
    send_requests = []
    recv_requests = []
    for iinterface in range(num_interfaces):
        npts = nibool_interfaces[iinterface]
        req = comm.Isend(
            # status, req = numba_mpi.isend(
            buffer_send_scalar[iinterface, :npts],
            dest=my_neighbors[iinterface],
            tag=11,
        )
        send_requests.append(req)
        # send_requests[iinterface] = req[0]
        req = comm.Irecv(
            # status, req = numba_mpi.irecv(
            buffer_recv_scalar[iinterface, :npts],
            source=my_neighbors[iinterface],
            tag=11,
        )
        recv_requests.append(req)
        # recv_requests[iinterface] = req[0]

    # wait for communications completion (recv)
    MPI.Request.Waitall(recv_requests)
    # numba_mpi.waitall(recv_requests)

    # adding contributions of neighbors, in the order of mpi_ranks
    sum_glob = np.zeros_like(array_glob)

    inds = np.argsort(my_neighbors)
    mask = my_neighbors[inds] < mpi_rank

    for i in inds[mask]:
        npts = nibool_interfaces[i]
        sum_glob[ibool_interfaces[i, :npts]] += buffer_recv_scalar[i, :npts]

    sum_glob[:] += array_glob[:]

    for i in inds[~mask]:
        npts = nibool_interfaces[i]
        sum_glob[ibool_interfaces[i, :npts]] += buffer_recv_scalar[i, :npts]

    # for i in range(num_interfaces):
    #     npts = nibool_interfaces[i]
    #     array_glob[ibool_interfaces[i, :npts]] += buffer_recv_scalar[i, :npts]

    # wait for communications completion (send)
    MPI.Request.Waitall(send_requests)
    # numba_mpi.waitall(send_requests)

    array_glob[:] = sum_glob[:]

    # return sum_glob


@numba.jit
def gll2glob(
    u_gll,
    nglob,
    ibool,  # [0:nspec,0:nzgll,0:nygll,0:nxgll]
):
    """
    u_gll[nspec,nzgll,nygll,nxgll]
    """
    nspec, ngllz, nglly, ngllx = ibool.shape
    u_glob = np.zeros(nglob, dtype=u_gll.dtype)
    # counts = np.zeros(nglob, dtype=ibool.dtype)

    for e in range(nspec):
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    # counts[idof] += 1
                    u_glob[idof] = u_glob[idof] + u_gll[e, k, j, i]
    # u_glob = u_glob / counts
    return u_glob


@numba.jit
def laplacian(
    u_glob,
    kappa,
    wgll,
    dlag_gll,
    ibool,
    dxsi_dx,
    dxsi_dy,
    dxsi_dz,
    deta_dx,
    deta_dy,
    deta_dz,
    dgam_dx,
    dgam_dy,
    dgam_dz,
    jacobian,
):
    """int(grad(phi_gll) * K * grad(u), dV)"""
    nspec, ngllz, nglly, ngllx = ibool.shape
    # assert u_glob.shape == out_glob.shape

    dtype = u_glob.dtype
    sl = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_xsil = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_etal = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_gaml = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    stif = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    out_glob = np.zeros(u_glob.shape, dtype=dtype)

    for e in range(nspec):

        # sl = u_glob[ibool[e, :, :, :]]
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    sl[k, j, i] = u_glob[idof]

        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Compute derivatives at elemental level
                    dsl_dxsi = 0
                    dsl_deta = 0
                    dsl_dgam = 0

                    for l in range(ngllx):
                        dsl_dxsi = dsl_dxsi + sl[k, j, l] * dlag_gll[l, i]
                        dsl_deta = dsl_deta + sl[k, l, i] * dlag_gll[l, j]
                        dsl_dgam = dsl_dgam + sl[l, j, i] * dlag_gll[l, k]

                    # Get derivatives informations
                    dxsi_dxl = dxsi_dx[e, k, j, i]
                    dxsi_dyl = dxsi_dy[e, k, j, i]
                    dxsi_dzl = dxsi_dz[e, k, j, i]
                    deta_dxl = deta_dx[e, k, j, i]
                    deta_dyl = deta_dy[e, k, j, i]
                    deta_dzl = deta_dz[e, k, j, i]
                    dgam_dxl = dgam_dx[e, k, j, i]
                    dgam_dyl = dgam_dy[e, k, j, i]
                    dgam_dzl = dgam_dz[e, k, j, i]
                    jacobianl = jacobian[e, k, j, i]

                    # Get physical derivatives
                    dsl_dxl = (
                        dsl_dxsi * dxsi_dxl + dsl_deta * deta_dxl + dsl_dgam * dgam_dxl
                    )
                    dsl_dyl = (
                        dsl_dxsi * dxsi_dyl + dsl_deta * deta_dyl + dsl_dgam * dgam_dyl
                    )
                    dsl_dzl = (
                        dsl_dxsi * dxsi_dzl + dsl_deta * deta_dzl + dsl_dgam * dgam_dzl
                    )

                    # Start product: Kij * du_dxi * dxsi_dxj
                    # kl = kappa_gll[e, k, j, i]
                    grad_xsil[k, j, i] = (
                        jacobianl
                        * kappa
                        * (dxsi_dxl * dsl_dxl + dxsi_dyl * dsl_dyl + dxsi_dzl * dsl_dzl)
                    )
                    grad_etal[k, j, i] = (
                        jacobianl
                        * kappa
                        * (deta_dxl * dsl_dxl + deta_dyl * dsl_dyl + deta_dzl * dsl_dzl)
                    )
                    grad_gaml[k, j, i] = (
                        jacobianl
                        * kappa
                        * (dgam_dxl * dsl_dxl + dgam_dyl * dsl_dyl + dgam_dzl * dsl_dzl)
                    )

        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Compute derivatives at elemental level
                    lapla_x = 0.0
                    lapla_y = 0.0
                    lapla_z = 0.0
                    for l in range(ngllx):
                        lapla_x = (
                            lapla_x + grad_xsil[k, j, l] * dlag_gll[i, l] * wgll[l]
                        )
                        lapla_y = (
                            lapla_y + grad_etal[k, l, i] * dlag_gll[j, l] * wgll[l]
                        )
                        lapla_z = (
                            lapla_z + grad_gaml[l, j, i] * dlag_gll[k, l] * wgll[l]
                        )

                    # Stiffness
                    stif[k, j, i] = (
                        wgll[j] * wgll[k] * lapla_x
                        + wgll[i] * wgll[k] * lapla_y
                        + wgll[i] * wgll[j] * lapla_z
                    )

        # Go back to dof
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    out_glob[idof] += stif[k, j, i]

    # numba_mpi.barrier()
    return out_glob


# ====== parameters
parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int)  # "misfit.h5"
parser.add_argument("mesh_dir")  # "search.txt"
parser.add_argument("model_dir")  # "search.pdf"
parser.add_argument("model_name")  # "search.pdf"
parser.add_argument("smooth_length", type=float)  # one sigma
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

dt = 1.0 / nt
smooth_length /= R_EARTH_KM  # non-dimensionalize as SEM
kappa = 0.5 * smooth_length**2

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

print(f'{mpi_rank=}, {mesh_mpi["my_neighbors"]=}')

comm.Barrier()
if mpi_rank == 0:
    print("====== finish reading mesh/model ", time.time() - tic)
    sys.stdout.flush()

dv_gll = sem_mesh_get_vol_gll(mesh)
dv_glob = gll2glob(
    dv_gll,
    nglob,
    ibool,
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
xgll, wgll = gll_library.zwgljd(NGLLX, GAUSSALPHA, GAUSSBETA)
# if mpi_rank == 0: print(f"{xgll=}, {wgll=}")

# Get derivative matrices
dlag_gll = np.zeros((NGLLX, NGLLX))  # , dtype=u_glob.dtype)
# dlag_wgll = np.zeros((NGLLX, NGLLX), dtype=u_glob.dtype)
for i in range(NGLLX):
    for j in range(NGLLX):
        # dLag_i/dx(xgll[j])
        dlag_gll[i, j] = gll_library.lagrange_deriv_gll(i, j, xgll, NGLLX)
        # dlag_wgll[i, j] = dlag_gll[i, j] * wgll[j]

GLL_Data = {"xgll": xgll, "wgll": wgll, "dlag_gll": dlag_gll}


# trapzoidal method
# M * (u_{n+1} - u_n) = dt * K * (u_{n+1} + u_n) / 2
# (M - dt/2 * K) * u_{n+1} = (M + dt/2 * K) * u_n
# (1 - dt/2 * M^{-1} * K) * u_{n+1} = (1 + dt/2 * M^{-1} * K) * u_n

def Kx(x_glob):
    kx_glob = -1.0 * laplacian(
        x_glob,
        kappa,
        GLL_Data["wgll"],
        GLL_Data["dlag_gll"],
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

# trapzoidal method
# M * (u_{n+1} - u_n) = dt * K * (u_{n+1} + u_n) / 2
# (M - dt/2 * K) * u_{n+1} = (M + dt/2 * K) * u_n
# (1 - dt/2 * M^{-1} * K) * u_{n+1} = (1 + dt/2 * M^{-1} * K) * u_n
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
    print("====== after smoothing ", time.time() - tic)
    # print(f"{du_glob=}, {du_glob.dtype=}")

out_file = "%s/proc%06d_reg1_%s.bin" % (out_dir, mpi_rank, model_name)
with FortranFile(out_file, "w") as f:
    # f.write_record(np.array(u_glob[ibool], dtype="f4"))
    f.write_record(np.array(u[ibool], dtype="f4"))