#!/bin/bash

# configure files: Par_file, CMTSOLUTION, setup/*.h.in

wkdir=$(pwd)

sem_config_dir=$wkdir/${1:-sem_config}
sem_build_dir=$wkdir/${2:-specfem3d_globe}

# prepare sem_build_dir
cd $sem_build_dir/
mkdir DATABASES_MPI OUTPUT_FILES

cd $sem_build_dir/DATA
ln -sf $sem_config_dir/DATA/Par_file .
ln -sf $sem_config_dir/DATA/CMTSOLUTION .

cd $sem_build_dir/setup
ln -sf $sem_config_dir/setup/*.h.in .

# build
cd $sem_build_dir

#./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-fc=ifort -O3"
./configure FC=mpif90 MPIFC=mpif90 FCFLAGS=" -O3"
make clean
make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk |\
  tee make.log