#!/bin/bash

# create SEM build directory and compile the program

sem_source_dir=${1:?[arg] need sem_source_dir}
sem_config_dir=${2:-sem_config}
sem_build_dir=${3:-specfem3d_globe}

sem_source_dir=$(readlink -f $sem_source_dir)
sem_config_dir=$(readlink -f $sem_config_dir)
sem_build_dir=$(readlink -f $sem_build_dir)

if [ ! -d "$sem_source_dir" ]
then
  echo "[ERROR] sem_source_dir does not exits!"
  exit 1
fi
if [ ! -d "$sem_config_dir" ]
then
  echo "[ERROR] sem_config_dir does not exits!"
  exit 1
fi

#====== prepare sem_build_dir
if [ -d "$sem_build_dir" ]
then
  echo "[WARN] sem_build_dir exits, delete!"
  rm -rf $sem_build_dir
fi
mkdir $sem_build_dir

cd $sem_build_dir
mkdir DATABASES_MPI OUTPUT_FILES

# link src/, config*
ln -s ${sem_source_dir}/src ./
ln -s ${sem_source_dir}/config* ./
ln -s ${sem_source_dir}/flags.guess ./
ln -s ${sem_source_dir}/install-sh ./
ln -s ${sem_source_dir}/Makefile.in ./

# link DATA/*
mkdir DATA; cd DATA
ln -s ${sem_source_dir}/DATA/* ./
cd ../

# copy setup directory
cp -r ${sem_source_dir}/setup ./

#====== use sem_config_dir
cd $sem_build_dir/DATA
ln -sf $sem_config_dir/DATA/Par_file .
ln -sf $sem_config_dir/DATA/CMTSOLUTION .

cd $sem_build_dir/setup
ln -sf $sem_config_dir/setup/*.h.in .

#====== build SEM 
cd $sem_build_dir

#./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-fc=ifort -O3"
./configure FC=mpif90 MPIFC=mpif90 FCFLAGS=" -O3"

make clean
make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk \
  > >(tee make.out) 2> >(tee make.err >&2)

#END