#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import warnings
import argparse
import numpy as np
from scipy.io import FortranFile

#====== user input
parser = argparse.ArgumentParser(description='Make vtk file from movie_volume files.')
parser.add_argument(dest='mesh_dir', type=str, help="directory containing movie3D_x|y|z|elements|DIZ00****.bin")
parser.add_argument(dest='iproc', type=int, help="iproc (< nproc)")
parser.add_argument(dest='it', type=int, help="iteration time step")
parser.add_argument(dest='wavefield_name', type=str, help="wavefield name, e.g. DIZ")
parser.add_argument(dest='out_dir', type=str, help="output directory")
args = parser.parse_args()
print(args)

mesh_dir = args.mesh_dir
iproc = args.iproc
it = args.it
wavefield_name = args.wavefield_name
out_dir = args.out_dir

#mesh_dir = "DATABASES_MPI/"
#iproc = 0
#it = 1200
#wavefield_name = 'DIZ'
vtk_file = "%s/proc%06d_movie3D_%s%06d.vtk"%(out_dir,iproc,wavefield_name,it)

# read in x-y-z coordinates
mesh = {}
for comp in ['x','y','z']:
    filename = "%s/proc%06d_movie3D_%s.bin"%(mesh_dir,iproc,comp)
    with FortranFile(filename, 'r') as f:
        mesh[comp] = f.read_ints(dtype='f4')

# read in cell topology: 8 corners of each hexahedron
filename = "%s/proc%06d_movie3D_elements.bin"%(mesh_dir,iproc)
with FortranFile(filename, 'r') as f:
    cells = []
    while True:
        try:
            cells.append(f.read_ints(dtype='i4'))
        except Exception as e:
            print(e)
            break

# read in wavefield snapshot
filename = "%s/proc%06d_movie3D_%s%06d.bin"%(mesh_dir,iproc,wavefield_name,it)
with FortranFile(filename, 'r') as f:
    with FortranFile(filename, 'r') as f:
        wavefield = f.read_ints(dtype='f4')

# write out VTK file
num_points = mesh['x'].size
num_cells = len(cells)

with open(vtk_file, 'w') as f:
  f.write("# vtk DataFile Version 2.0\n")
  f.write("movie_volume\n")
  f.write("ASCII\n")
  f.write("DATASET UNSTRUCTURED_GRID\n")
  f.write("POINTS %d float\n"%(num_points))
  for ip in range(num_points): 
    f.write("%+15.7E %+15.7E %+15.7E\n"%(mesh['x'][ip],mesh['y'][ip],mesh['z'][ip]))
  f.write("CELLS %d %d\n"%(num_cells, num_cells+8*num_cells)) # HEXAHEDRON
  for ie in range(num_cells):
    f.write("8 %d %d %d %d %d %d %d %d\n"%(cells[ie][0],cells[ie][1],cells[ie][2],cells[ie][3],cells[ie][4],cells[ie][5],cells[ie][6],cells[ie][7]))
  f.write("CELL_TYPES %d\n"%(num_cells)) # HEXAHEDRON
  for ie in range(num_cells):
    f.write("12\n")
  f.write("POINT_DATA %d\n"%(num_points))
  f.write("SCALARS %s float 1\n"%(wavefield_name))
  f.write("LOOKUP_TABLE default\n")
  for ip in range(num_points): 
    f.write("%+15.7E\n"%(wavefield[ip]))
