#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""create horizontal slice of SEM model at a given depth"""
import sys
import time
import argparse

import numpy as np
from scipy.io import FortranFile
import numba

from meshfem3d_utils import sem_mesh_read, R_EARTH_KM

# ======
parser = argparse.ArgumentParser()

parser.add_argument("--nproc", help="number of slices", type=int)
parser.add_argument("--mesh_dir", help="input mesh directory of proc*_solver_data.bin")
parser.add_argument("--model_dir", help="input model directory of proc*_reg1_[model_tag].bin")
parser.add_argument("--model_name", help="tag for velocity model, e.g. vsv")
parser.add_argument("--dt", help="time step in second, e.g. 0.1", type=float)

args = parser.parse_args()
print(args)

nproc = args.nproc  
mesh_dir = args.mesh_dir
model_dir = args.model_dir
model_name = args.model_name
dt = args.dt

# nproc = int(sys.argv[1])
# mesh_dir = str(sys.argv[2])  # <mesh_dir>/proc******_external_mesh.bin
# model_dir = str(sys.argv[3])
# model_name = str(sys.argv[4])
# dt = float(sys.argv[5])


@numba.jit(nopython=True)
def _calc_CFL(ibool, xyz_glob, vel, nspec, dt):
    CFL = np.zeros(nspec)
    for ispec in range(nspec):
        iglob = ibool[ispec, :, :, :].flatten()
        xyz_gll = xyz_glob[iglob, :]

        # min_dist_gll = np.inf
        # n = iglob.size
        # for i in range(n):
        #     for j in range(i + 1, n):
        #         dist = np.sum((xyz_gll[i, :] - xyz_gll[j, :])**2)
        #         if dist < min_dist_gll:
        #             min_dist_gll = dist
        # min_dist_gll = min_dist_gll**0.5 * R_EARTH_KM

        dist2 = (
            np.sum((xyz_gll[:, None, :] - xyz_gll[None, :, :]) ** 2, axis=2) ** 0.5
            * R_EARTH_KM
        )
        np.fill_diagonal(dist2, np.inf)
        min_dist_gll = np.min(dist2)
        
        CFL[ispec] = dt * np.max(vel[ispec, :, :, :]) / min_dist_gll
    return CFL


for iproc in range(nproc):

    mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
    mesh_data = sem_mesh_read(mesh_file)
    nspec = mesh_data["nspec"]
    ibool = mesh_data["ibool"]
    xyz_glob = mesh_data["xyz_glob"]
    gll_dims = mesh_data["gll_dims"]

    model_file = "%s/proc%06d_reg1_%s.bin" % (model_dir, iproc, model_name)
    with FortranFile(model_file, "r") as f:
        vel = np.reshape(f.read_reals(dtype="f4"), gll_dims)

    # tic = time.time()
    CFL = _calc_CFL(ibool, xyz_glob, vel, nspec, dt)
    # toc = time.time()
    # print("iproc=%d, time=%f" % (iproc, toc-tic))

    print(
        "[proc%03d] min/max vel= %f %f, min/max CFL= %f %f"
        % (iproc, np.min(vel), np.max(vel), np.min(CFL), np.max(CFL))
    )
