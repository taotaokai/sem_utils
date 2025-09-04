#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import os
import argparse
import numpy as np
from scipy.io import FortranFile
import pyvista as pv

from meshfem3d_utils import sem_mesh_read, NGLLX, NGLLY, NGLLZ

#====== user input
parser = argparse.ArgumentParser(description='Make sem mesh vtk file')
parser.add_argument(dest='procnum_begin', type=int, help="begin number of iproc (>= 0)")
parser.add_argument(dest='procnum_end', type=int, help="end number of iproc (< nproc)")
parser.add_argument(dest='mesh_dir', type=str, help="DATABSAES_MPI/")
parser.add_argument(dest='model_dir', type=str, help="DATABSAES_MPI/")
parser.add_argument(dest='model_tag', type=str, help="vsv")
parser.add_argument(dest='vtk_file', type=str, help="output .vtk file name")
args = parser.parse_args()
print(args)

#======
point_xyz_total = np.zeros((0,3), dtype=np.float32)
cell_topo_total = np.zeros((0,8), dtype=np.uint32)
point_data_total = np.zeros((0), dtype=np.float32)
point_num_total = 0
cell_num_total = 0

# indices of the HEX points in GLL element
corner_inds = [
         0  + NGLLX * (      0 + NGLLY * (      0)),  
   NGLLX-1  + NGLLX * (      0 + NGLLY * (      0)),  
   NGLLX-1  + NGLLX * (NGLLY-1 + NGLLY * (      0)),  
         0  + NGLLX * (NGLLY-1 + NGLLY * (      0)),  
         0  + NGLLX * (      0 + NGLLY * (NGLLZ-1)),  
   NGLLX-1  + NGLLX * (      0 + NGLLY * (NGLLZ-1)),  
   NGLLX-1  + NGLLX * (NGLLY-1 + NGLLY * (NGLLZ-1)),  
         0  + NGLLX * (NGLLY-1 + NGLLY * (NGLLZ-1)),  
]

for iproc in range(args.procnum_begin, args.procnum_end + 1):

    print('# iproc = ', iproc)

    prname = f'proc{iproc:06d}_reg1'
    mesh_file = os.path.join(args.mesh_dir, f'{prname}_solver_data.bin')
    mesh_data = sem_mesh_read(mesh_file)

    gll_dims = mesh_data["gll_dims"]
    model_file = os.path.join(args.model_dir, f'{prname}_{args.model_tag}.bin')
    with FortranFile(model_file, "r") as f:
        # here gll_dims = [NSPEC, NGLLZ, NGLLY, NGLLX] which is C convention of
        # Fortran array of [NGLLX, NGLLY, NGLLZ, NSPEC]
        data = f.read_ints(dtype="f4")
        model_gll = np.reshape(data, gll_dims)

    # vtk point/cell data
    nspec = mesh_data['nspec']
    nglob = mesh_data['nglob']
    ibool = mesh_data['ibool']
    idoubling = mesh_data['idoubling']
    xyz_glob = mesh_data['xyz_glob']
    # convert [NSPEC, NGLLZ,NGLLY,NGLLX] to [NSPEC, NGLLZ*NGLLY*NGLLX]
    ibool = ibool.reshape((nspec, -1))
    model_gll = model_gll.reshape((nspec, -1))

    print('nglob,nspec = ', nglob, nspec)
    new_iglob = np.ones(nglob, dtype=bool)
    num_iglob = np.zeros(nglob, dtype=np.uint32)
    point_xyz = np.zeros((nglob,3), dtype=np.float32)
    point_data = np.zeros(nglob, dtype=np.float32)
    cell_topo = np.zeros((nspec,8), dtype=np.uint32)  # 8 corners of each cell

    point_num = 0 # point count
    cell_num = 0 # cell count
    for ispec in range(nspec):
        ipoint_cell = 0 # point count in one cell
        # 8 corners of each hex cell
        iglob_list = ibool[ispec, corner_inds]
        for i, iglob in enumerate(iglob_list):
            if new_iglob[iglob]: # register new iglob into point_data
              new_iglob[iglob] = False
              num_iglob[iglob] = point_num + point_num_total
              point_xyz[point_num, :] = xyz_glob[iglob, :]
              # NOTE point values shared by neighouring element will be ignored 
              point_data[point_num] = model_gll[ispec, corner_inds[i]]
              point_num = point_num + 1
            cell_topo[cell_num, ipoint_cell] = num_iglob[iglob]
            ipoint_cell = ipoint_cell + 1
        cell_num = cell_num + 1

    point_num_total = point_num_total + point_num
    cell_num_total = cell_num_total + cell_num

    point_xyz_total = np.vstack((point_xyz_total, point_xyz[0:point_num, :]))
    cell_topo_total = np.vstack((cell_topo_total, cell_topo[0:cell_num, :]))
    point_data_total = np.hstack((point_data_total, point_data[0:point_num]))

# write out vtk
grid = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: cell_topo_total}, point_xyz_total)
grid.point_data[args.model_tag] = point_data_total
grid.save(args.vtk_file)
