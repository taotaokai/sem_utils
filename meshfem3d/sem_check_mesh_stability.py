#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""create horizontal slice of SEM model at a given depth"""
import sys

import numpy as np
from scipy.io import FortranFile

from meshfem3d_utils import sem_mesh_read, R_EARTH_KM

# ======
nproc = int(sys.argv[1])
mesh_dir = str(sys.argv[2])  # <mesh_dir>/proc******_external_mesh.bin
model_dir = str(sys.argv[3])
model_name = str(sys.argv[4])
dt = float(sys.argv[5])

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

    # --- determine element size (approximately)
    CFL = np.zeros(nspec)
    for ispec in range(nspec):
        iglob = ibool[ispec, :, :, :].reshape(-1)
        xyz_gll = xyz_glob[iglob, :]
        dist2 = (
            np.sum((xyz_gll[:, None, :] - xyz_gll[None, :, :]) ** 2, axis=0) ** 0.5
            * R_EARTH_KM
        )
        np.fill_diagonal(dist2, np.nan)
        min_dist_gll = np.nanmin(dist2)
        CFL[ispec] = dt * np.max(vel[ispec, :, :, :]) / min_dist_gll

    print(
        "[proc%03d] min/max vel= %f %f, min/max CFL= %f %f"
        % (iproc, np.min(vel), np.max(vel), np.min(CFL), np.max(CFL))
    )
