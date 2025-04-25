#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import warnings
import numpy as np
import argparse
from meshfem3d_utils import  sem_mesh_read,sem_boundary_mesh_read,NGLLX,NGLLY,NGLLZ

#====== user input
parser = argparse.ArgumentParser(description='Make vtk file from the teleseismic_boundary.bin files.')
parser.add_argument(dest='procnum_begin', type=int, help="begin number of iproc (>= 0)")
parser.add_argument(dest='procnum_end', type=int, help="end number of iproc (< nproc)")
parser.add_argument(dest='mesh_dir', type=str, help="directory containing solver_data.bin and teleseismic_boundary.bin")
parser.add_argument('--output', type=argparse.FileType('w', encoding='UTF-8'), default='teleseismic_boundary.vtk', dest='vtk_file')
args = parser.parse_args()
print(args)

#procnum_begin = 0
#procnum_end = 4
#mesh_dir = 'mesh/DATABASES_MPI/'
#vtk_file = 'teleseismic_boundary.vtk'

#====== 
point_data_total = np.zeros((3,0), dtype=np.float)
cell_data_total = np.zeros((4,0), dtype=np.uint32)
point_num_total = 0
cell_num_total = 0
 
for iproc in range(args.procnum_begin, args.procnum_end + 1):

  print('# iproc = ', iproc)

  prname = args.mesh_dir + '/proc%06d_reg1_'%(iproc)
  mesh_file = prname + 'solver_data.bin'
  mesh_data = sem_mesh_read(mesh_file)
  boundary_file = prname + 'teleseismic_boundary.bin'
  boundary_data = sem_boundary_mesh_read(boundary_file)

  print(mesh_file)
  
  # vtk point/cell data
  nspec = mesh_data['nspec']
  nglob = mesh_data['nglob']
  ibool = mesh_data['ibool']
  xyz_glob = mesh_data['xyz_glob']

  print('nglob,nspec = ', nglob, nspec)
  new_iglob = np.ones(nglob, dtype=np.bool)
  num_iglob = np.zeros(nglob, dtype=np.uint32)
  point_data = np.zeros((3,nglob), dtype=np.float)
  cell_data = np.zeros((4,nspec), dtype=np.uint32)

  point_num = 0 # point count
  cell_num = 0 # cell count
  for position in ['xmin','xmax','ymin','ymax','zmin']:
    tag = 'ibelm_teleseismic_' + position 
    for ispec in boundary_data[tag]:
      ipoint_cell = 0 # point count in one cell
      if position == 'xmin': iglob_list = [ibool[0,0,0,ispec-1], ibool[0,0,NGLLZ-1,ispec-1],ibool[0,NGLLY-1,NGLLZ-1,ispec-1],ibool[0,NGLLY-1,0,ispec-1]]
      if position == 'xmax': iglob_list = [ibool[NGLLX-1,0,0,ispec-1], ibool[NGLLX-1,NGLLY-1,0,ispec-1],ibool[NGLLX-1,NGLLY-1,NGLLZ-1,ispec-1],ibool[NGLLX-1,0,NGLLZ-1,ispec-1]]
      if position == 'ymin': iglob_list = [ibool[0,0,0,ispec-1], ibool[NGLLX-1,0,0,ispec-1],ibool[NGLLX-1,0,NGLLZ-1,ispec-1],ibool[0,0,NGLLZ-1,ispec-1]]
      if position == 'ymax': iglob_list = [ibool[0,NGLLY-1,0,ispec-1], ibool[NGLLX-1,NGLLY-1,0,ispec-1],ibool[NGLLX-1,NGLLY-1,NGLLZ-1,ispec-1],ibool[0,NGLLY-1,NGLLZ-1,ispec-1]]
      if position == 'zmin': iglob_list = [ibool[0,0,0,ispec-1], ibool[NGLLX-1,0,0,ispec-1],ibool[NGLLX-1,NGLLY-1,0,ispec-1],ibool[0,NGLLY-1,0,ispec-1]]
      iglob_list = [ i-1 for i in iglob_list]
      for iglob in iglob_list:
        if new_iglob[iglob]: # register unseen iglob into point_data
          new_iglob[iglob] = False
          num_iglob[iglob] = point_num + point_num_total
          point_data[:,point_num] = xyz_glob[:,iglob]
          point_num = point_num + 1
        cell_data[ipoint_cell,cell_num] = num_iglob[iglob]
        ipoint_cell = ipoint_cell + 1
      cell_num = cell_num + 1

  print('point_num,cell_num =', point_num, cell_num)

  point_num_total = point_num_total + point_num
  cell_num_total = cell_num_total + cell_num

  point_data_total = np.hstack((point_data_total, point_data[:,0:point_num]))
  cell_data_total = np.hstack((cell_data_total, cell_data[:,0:cell_num]))

# write out vtk
with args.vtk_file as f:
  f.write("# vtk DataFile Version 2.0\n")
  f.write("boundary mesh\n")
  f.write("ASCII\n")
  f.write("DATASET UNSTRUCTURED_GRID\n")
  f.write("POINTS %d float\n"%(point_num_total))
  for ip in range(point_num_total): 
      f.write("%+15.7E %+15.7E %+15.7E\n"%(point_data_total[0,ip],point_data_total[1,ip],point_data_total[2,ip]))
  f.write("CELLS %d %d\n"%(cell_num_total, cell_num_total+4*cell_num_total)) # QUAD
  for ie in range(cell_num_total):
    f.write("4 %d %d %d %d\n"%(cell_data_total[0,ie],cell_data_total[1,ie],cell_data_total[2,ie],cell_data_total[3,ie]))
  f.write("CELL_TYPES %d\n"%(cell_num_total)) # QUAD
  for ie in range(cell_num_total):
    f.write("9\n")
