#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys
import argparse
import numpy as np
from scipy.io import FortranFile
import pyvista as pv

from meshfem3d_utils import sem_mesh_read, NGLLX, NGLLY, NGLLZ

#====== user input
parser = argparse.ArgumentParser(description='Make vtk file from the teleseismic_boundary.bin files.')
parser.add_argument(dest='procnum_begin', type=int, help="begin number of iproc (>= 0)")
parser.add_argument(dest='procnum_end', type=int, help="end number of iproc (< nproc)")
parser.add_argument(dest='mesh_dir', type=str, help="directory containing solver_data.bin and teleseismic_boundary.bin")
parser.add_argument(dest='vtk_file', type=str, help="output .vtk file name")
args = parser.parse_args()
print(args)

#======
point_xyz_total = np.zeros((0,3), dtype=np.float32)
cell_topo_total = np.zeros((0,8), dtype=np.uint32)
cell_data_total = np.zeros((0), dtype=np.int32)
point_num_total = 0
cell_num_total = 0

for iproc in range(args.procnum_begin, args.procnum_end + 1):

    print('# iproc = ', iproc)

    prname = args.mesh_dir + '/proc%06d_reg1_'%(iproc)
    mesh_file = prname + 'solver_data.bin'
    mesh_data = sem_mesh_read(mesh_file)

    print(mesh_file)

    # vtk point/cell data
    nspec = mesh_data['nspec']
    nglob = mesh_data['nglob']
    ibool = mesh_data['ibool']
    idoubling = mesh_data['idoubling']
    xyz_glob = mesh_data['xyz_glob']

    print('nglob,nspec = ', nglob, nspec)
    new_iglob = np.ones(nglob, dtype=bool)
    num_iglob = np.zeros(nglob, dtype=np.uint32)
    point_xyz = np.zeros((nglob,3), dtype=np.float32)
    cell_topo = np.zeros((nspec,8), dtype=np.uint32)  # 8 corners of each cell

    point_num = 0 # point count
    cell_num = 0 # cell count
    for ispec in range(nspec):
        ipoint_cell = 0 # point count in one cell
        iglob_list = [  # 8 corners of each hex cell
            ibool[ispec,       0,       0, 0      ],
            ibool[ispec,       0,       0, NGLLX-1],
            ibool[ispec,       0, NGLLY-1, NGLLX-1],
            ibool[ispec,       0, NGLLY-1,       0],
            ibool[ispec, NGLLZ-1,       0, 0      ],
            ibool[ispec, NGLLZ-1,       0, NGLLX-1],
            ibool[ispec, NGLLZ-1, NGLLY-1, NGLLX-1],
            ibool[ispec, NGLLZ-1, NGLLY-1,       0],
        ]
        for iglob in iglob_list:
            if new_iglob[iglob]: # register unseen iglob into point_data
              new_iglob[iglob] = False
              num_iglob[iglob] = point_num + point_num_total
              point_xyz[point_num, :] = xyz_glob[iglob, :]
              point_num = point_num + 1
            cell_topo[cell_num, ipoint_cell] = num_iglob[iglob]
            ipoint_cell = ipoint_cell + 1
        cell_num = cell_num + 1

    point_num_total = point_num_total + point_num
    cell_num_total = cell_num_total + cell_num

    point_xyz_total = np.vstack((point_xyz_total, point_xyz[0:point_num, :]))
    cell_topo_total = np.vstack((cell_topo_total, cell_topo[0:cell_num, :]))
    cell_data_total = np.hstack((cell_data_total, np.arange(nspec)))

# write out vtk
grid = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: cell_topo_total}, point_xyz_total)
grid.cell_data['ispec'] = cell_data_total
grid.save(args.vtk_file)
