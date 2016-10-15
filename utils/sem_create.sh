#!/bin/bash

# create SEM build directory and compile the program

sem_source_dir=${1:?[arg] need sem_source_dir}
sem_setup_dir=${2:?[arg] need setup dir(for constants.h.in)}
sem_data_dir=${3:?[arg] need DATA dir(for DATA/Par_file,STATIONS,CMTSOLUTION)}
sem_build_dir=${4:?[arg] need build dir(e.g. specfem3d_globe)}

sem_source_dir=$(readlink -f $sem_source_dir)
sem_setup_dir=$(readlink -f $sem_setup_dir)
sem_data_dir=$(readlink -f $sem_data_dir)
sem_build_dir=$(readlink -f $sem_build_dir)

if [ ! -d "$sem_source_dir" ]
then
  echo "[ERROR] sem_source_dir does not exits!"
  exit 1
fi
if [ ! -d "$sem_setup_dir" ]
then
  echo "[ERROR] sem_setup_dir does not exits!"
  exit 1
fi
if [ ! -d "$sem_data_dir" ]
then
  echo "[ERROR] sem_data_dir does not exits!"
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
ln -s ${sem_source_dir}/config.guess ./
ln -s ${sem_source_dir}/config.sub ./
ln -s ${sem_source_dir}/configure ./
ln -s ${sem_source_dir}/configure.ac ./
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
rm -f Par_file CMTSOLUTION
cp $sem_data_dir/Par_file .
cp $sem_data_dir/CMTSOLUTION .

cd $sem_build_dir/setup
rm -f constants.h.in
cp $sem_setup_dir/constants.h.in .

#====== build SEM 
cd $sem_build_dir

#./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-fc=ifort -O3"
./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-O3" CC=mpicc CFLAGS="-O3"

make clean
make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk \
  > >(tee make.out) 2> >(tee make.err >&2)

#END