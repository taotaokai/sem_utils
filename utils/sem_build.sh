#!/bin/bash

# configure files: Par_file, CMTSOLUTION, setup/*.h.in

wkdir=$(pwd)

config_dir=$wkdir/${1:-sem_config}
build_dir=$wkdir/${2:-specfem3d_globe}

# prepare build_dir
cd $build_dir/
mkdir DATABASES_MPI OUTPUT_FILES

cd $build_dir/DATA
ln -sf $config_dir/DATA/Par_file .
ln -sf $config_dir/DATA/CMTSOLUTION .

cd $build_dir/setup
ln -sf $config_dir/setup/*.h.in .

# build
cd $build_dir

./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-O3"
make clean
make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk |\
  tee make.log