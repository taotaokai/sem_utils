#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""create horizontal slice of SEM model at a given depth"""
import sys
import warnings
import time
import argparse

import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import interpn
from mpi4py import MPI
import numba

# from meshfem3d_constants import NGLLX, NGLLY, NGLLZ, GAUSSALPHA, GAUSSBETA
from meshfem3d_constants import (
    IFLAG_CRUST,
    IFLAG_80_MOHO,
    IFLAG_220_80,
    IFLAG_670_220,
    IFLAG_DUMMY,
)

from meshfem3d_utils import (
    get_gll_weights,
    lagrange_poly,
    interp1d_linear,
    interp_model_gll,
    sem_mesh_read,
    sem_mesh_locate_points,
)

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
    "--method", default='linear', choices=['gll', 'linear'], help="Choose interpolation method (gll, linear)"
)

args = parser.parse_args()
print(args)

nproc_source = args.nproc_source
mesh_dir_source = args.mesh_dir_source
model_dir_source = args.model_dir_source

nproc_target = args.nproc_target
mesh_dir_target = args.mesh_dir_target
model_dir_target = args.model_dir_target

model_tags = args.model_tags
nmodel = len(model_tags)

# merge regions
idoubling_merge = []
# In SETibet case, since I use a velocity gradient across Moho and no mesh boundary at Moho, treat IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80 as the same region
# idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220]
# idoubling_merge = [IFLAG_CRUST, IFLAG_80_MOHO, IFLAG_220_80]

# ====== interpolate
comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

# GLL nodes
zgll, wgll, dlag_dzgll = get_gll_weights()

# --- loop over each slice of target SEM mesh
# for iproc_target in range(nproc_target):
for iproc_target in range(mpi_rank, nproc_target, mpi_size):
    # print("====== iproc_target ", iproc_target)
    # sys.stdout.flush()

    # read in target SEM mesh
    mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir_target, iproc_target)
    mesh_data_target = sem_mesh_read(mesh_file)
    nspec_target = mesh_data_target["nspec"]
    ibool_target = mesh_data_target["ibool"]
    idoubling_target = mesh_data_target["idoubling"]
    xyz_glob_target = mesh_data_target["xyz_glob"]
    # print(xyz_glob_target.shape)

    # merge regions if required
    idx_merge = np.zeros(nspec_target, dtype="bool")
    for ii in idoubling_merge:
        idx_merge = idx_merge | (idoubling_target == ii)
    idoubling_target[idx_merge] = IFLAG_DUMMY

    # xyz points to locate
    # xyz_target = xyz_glob_target[ibool_target.ravle(), :]
    # idoubling_ext = np.zeros(ibool_target.shape, dtype="int") + idoubling_target
    # idoubling_ext = idoubling_ext.ravel()
    xyz_target = xyz_glob_target
    idoubling_ext = -1

    # array of final results
    npoints = xyz_target.shape[0]
    status_target = np.zeros(npoints, dtype="int")
    status_target[:] = -1
    misloc_target = np.zeros(npoints)
    misloc_target[:] = np.inf
    misratio_target = np.zeros(npoints)
    model_target = np.zeros((nmodel, npoints))

    # -- loop over each slice of source SEM mesh
    for iproc_source in range(nproc_source):
        tic = time.time()

        # read in source SEM mesh
        mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir_source, iproc_source)
        mesh_data_source = sem_mesh_read(mesh_file)

        # merge regions if required
        idoubling_source = mesh_data_source["idoubling"]
        idx_merge = np.zeros(mesh_data_source["nspec"], dtype="bool")
        for ii in idoubling_merge:
            idx_merge = idx_merge | (idoubling_source == ii)
        idoubling_source[idx_merge] = IFLAG_DUMMY

        # read in source model
        gll_dims = mesh_data_source["gll_dims"]
        source_model_gll = np.zeros((nmodel,) + gll_dims)
        for imodel in range(nmodel):
            model_tag = model_tags[imodel]
            model_file = "%s/proc%06d_reg1_%s.bin" % (
                model_dir_source,
                iproc_source,
                model_tag,
            )
            with FortranFile(model_file, "r") as f:
                # here gll_dims = [NSPEC, NGLLZ, NGLLY, NGLLX] which is C convention of
                # Fortran array of [NGLLX, NGLLY, NGLLZ, NSPEC]
                source_model_gll[imodel, :, :, :, :] = np.reshape(
                    f.read_ints(dtype="f4"), gll_dims
                )

        # locate target points
        status_all, ispec_all, uvw_all, misloc_all, misratio_all = (
            sem_mesh_locate_points(
                mesh_data_source,
                xyz_target,
                idoubling_ext,
                kdtree_num_element=2.0,
                max_dist_ratio=2.0,
            )
        )

        # merge interpolation results of mesh slice (iproc_souce) into
        # the final results based on misloc and status

        # index selection for merge:
        # (not located inside an element yet) and
        # (located for the current mesh slice) and
        # ( smaller misloc or located inside an element in this mesh slice )
        ii = (
            (status_target != 1)
            & (status_all != -1)
            & ((misloc_all < misloc_target) | (status_all == 1))
        )

        status_target[ii] = status_all[ii]
        misloc_target[ii] = misloc_all[ii]
        misratio_target[ii] = misratio_all[ii]

        # NOTE avoid too many for loops reduces computation time
        ##for ipoint in range(npoints):
        # slower than index slicing but use less memory

        ipoint_select = np.nonzero(ii)[0]
        interp_model_gll(ipoint_select, zgll, ispec_all, uvw_all, model_gll, model_target, method=args.method)
        # for ipoint in ipoint_select:
        #     # interpolation weights
        #     if args.method == 'linear':
        #         wx = interp1d_linear(zgll, uvw_all[ipoint, 0])
        #         wy = interp1d_linear(zgll, uvw_all[ipoint, 1])
        #         wz = interp1d_linear(zgll, uvw_all[ipoint, 2])
        #     elif args.method == 'gll':
        #         wx = lagrange_poly(zgll, uvw_all[ipoint, 0])
        #         wy = lagrange_poly(zgll, uvw_all[ipoint, 1])
        #         wz = lagrange_poly(zgll, uvw_all[ipoint, 2])
        #     else:
        #         raise ValueError(f"Unknown interpolation method: {args.method}")
        #     # get interpolated values
        #     model_target[:, ipoint] = np.sum(
        #         source_model_gll[:, ispec_all[ipoint], :, :, :]
        #         * wx[None, None, None, :]
        #         * wy[None, None, :, None]
        #         * wz[None, :, None, None],
        #         axis=(1, 2, 3),
        #     )

        elapsed_time = time.time() - tic
        print(f"{iproc_target=:03d}, {iproc_source=:03d}, {elapsed_time=:8.3f} seconds")
        sys.stdout.flush()

    # rehape results
    gll_dims = mesh_data_target["gll_dims"]
    # status_gll_target = np.reshape(status_gll_target, gll_dims)
    # misloc_gll_target = np.reshape(misloc_gll_target, gll_dims)
    # misratio_gll_target = np.reshape(misratio_gll_target, gll_dims)
    # gll_dims = (nmodel,) + gll_dims
    # model_gll_target = np.reshape(model_gll_target, gll_dims)

    # save interpolated model
    for imodel in range(nmodel):
        model_tag = model_tags[imodel]
        model_file = "%s/proc%06d_reg1_%s.bin" % (
            model_dir_target,
            iproc_target,
            model_tag,
        )
        model_gll = model_target[imodel, ibool_target]
        with FortranFile(model_file, "w") as f:
            f.write_record(
                np.array(
                    np.ravel(model_gll),
                    dtype="f4",
                )
            )

    # save misloc, status
    model_file = "%s/proc%06d_reg1_status.bin" % (model_dir_target, iproc_target)
    with FortranFile(model_file, "w") as f:
        f.write_record(
            np.array(np.ravel(status_target[ibool_target], order="F"), dtype="f4")
        )

    model_file = "%s/proc%06d_reg1_misloc.bin" % (model_dir_target, iproc_target)
    with FortranFile(model_file, "w") as f:
        f.write_record(
            np.array(np.ravel(misloc_target[ibool_target], order="F"), dtype="f4")
        )

    model_file = "%s/proc%06d_reg1_misratio.bin" % (model_dir_target, iproc_target)
    with FortranFile(model_file, "w") as f:
        f.write_record(
            np.array(np.ravel(misratio_target[ibool_target], order="F"), dtype="f4")
        )