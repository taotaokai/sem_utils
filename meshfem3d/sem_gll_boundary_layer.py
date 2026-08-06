#!/usr/bin/env python
"""create boundary layer damping coefficient of SEM model"""
import argparse
import os

import numpy as np
from meshfem3d_utils import (
    read_sem_parfile,
    sem_mesh_read,
    sem_vec2xieta,
    write_gll_file,
)

# ====== user input
parser = argparse.ArgumentParser(description="Make sem gll files for boundary layer damping")
parser.add_argument("nproc", type=int, help="total number of slices")
parser.add_argument(
    "mesh_dir", type=str, help="{mesh_dir}/proc******_{region}_solver_data.bin"
)
parser.add_argument(
    "par_file", type=str, help="e.g. DATA/Par_file"
)
parser.add_argument("out_dir", type=str, help="output directory of gll files")
parser.add_argument(
    "--out_tag",
    default="taper",
    type=str,
    help="e.g. {out_dir}/proc****_{region}_{out_tag}.bin",
)
parser.add_argument(
    "--thickness",
    default=1.0,
    type=float,
    help="boundary layer thickness in degree",
)
parser.add_argument(
    "--taper",
    default=0.1,
    type=float,
    help="mininum in taper",
)
parser.add_argument(
    "--region",
    default="reg1",
    type=str,
    help="e.g. {mesh_dir}/proc****_{region}_solver_data.bin",
)
args = parser.parse_args()
print(args)

sem_params = read_sem_parfile(args.par_file)
mesh_central_lat = sem_params['CENTER_LATITUDE_IN_DEGREES']
mesh_central_lon = sem_params['CENTER_LONGITUDE_IN_DEGREES']
mesh_width_xi = sem_params['ANGULAR_WIDTH_XI_IN_DEGREES']
mesh_width_eta = sem_params['ANGULAR_WIDTH_ETA_IN_DEGREES']
mesh_gamma_rot = sem_params['GAMMA_ROTATION_AZIMUTH']
 

for iproc in range(args.nproc):

    print("# iproc = ", iproc)

    prname = f"proc{iproc:06d}_{args.region}"
    mesh_file = os.path.join(args.mesh_dir, f"{prname}_solver_data.bin")
    mesh_data = sem_mesh_read(mesh_file)

    gll_dims = mesh_data["gll_dims"]

    # vtk point/cell data
    nspec = mesh_data["nspec"]
    nglob = mesh_data["nglob"]
    ibool = mesh_data["ibool"]
    # idoubling = mesh_data["idoubling"]
    xyz_glob = mesh_data["xyz_glob"]

    xi, eta = sem_vec2xieta(mesh_central_lat, mesh_central_lon, mesh_gamma_rot, xyz_glob)

    # angular distance from edge  
    dxi = 0.5*mesh_width_xi - np.abs(xi)
    deta = 0.5*mesh_width_eta - np.abs(eta)
    
    dist2edge = np.zeros(xyz_glob.shape[:-1])

    # points inside the mesh
    ind = (dxi >= 0) & (deta >= 0)
    dist2edge[ind] = np.minimum(dxi[ind], deta[ind])

    taper = np.ones_like(dist2edge)
    ind = dist2edge <= args.thickness
    taper[ind] = args.taper + (1 - args.taper) * dist2edge[ind] / args.thickness

    # # points outside the mesh (as negative distance)
    # ind = (dxi < 0) & (deta > 0)
    # min_dist[ind] = dxi[ind]
    # ind = (dxi > 0) & (deta < 0)
    # min_dist[ind] = deta[ind]
    # ind = (dxi < 0) & (deta < 0)
    # min_dist[ind] = np.minimum(dxi[ind], deta[ind])
    # # multiply by radius to get distance
    # min_dist = min_dist * r[:, None]

    write_gll_file(args.out_dir, args.out_tag, iproc, taper[ibool], region_code=args.region, overwrite=True)