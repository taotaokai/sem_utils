#!/bin/bash

build_dir=specfem3d_globe
cd $build_dir 

make clean

./configure FC=mpif90 MPIFC=mpif90

make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk