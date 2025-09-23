#!/bin/bash

# create SEM build directory and compile the program

sem_source_dir=${1:?[arg] need sem_source_dir}
sem_setup_dir=${2:?[arg] need setup dir(for constants.h.in)}
sem_data_dir=${3:?[arg] need DATA dir(for DATA/Par_file,STATIONS,CMTSOLUTION)}
sem_build_dir=${4:?[arg] need build dir(e.g. specfem3d_globe)}
processor_type=${5:-cpu}

#mpif90=${5:?[arg]need mpif90 executable}
#mpicc=${6:?[arg]need mpicc executable}

sem_source_dir=$(readlink -f $sem_source_dir)
sem_setup_dir=$(readlink -f $sem_setup_dir)
sem_data_dir=$(readlink -f $sem_data_dir)
sem_build_dir=$(readlink -f $sem_build_dir)

if [ x${sem_source_dir} = x${sem_build_dir} ]
then
  echo "[ERROR] sem_source_dir cannot be the same as sem_build_dir!"
  exit -1
fi
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
if [ ! -f "$sem_data_dir/CMTSOLUTION" ]
then
  echo "[ERROR] sem_data_dir/DATA/CMTSOLUTION does not exits!"
  exit 1
fi
if [ ! -f "$sem_data_dir/STATIONS" ]
then
  echo "[ERROR] sem_data_dir/DATA/STATIONS does not exits!"
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
mkdir DATA DATABASES_MPI OUTPUT_FILES

# link src/, config*
ln -s ${sem_source_dir}/src ./
ln -s ${sem_source_dir}/config.guess ./
ln -s ${sem_source_dir}/config.sub ./
ln -s ${sem_source_dir}/configure ./
ln -s ${sem_source_dir}/configure.ac ./
ln -s ${sem_source_dir}/flags.guess ./
ln -s ${sem_source_dir}/Makefile.in ./
# link DATA/*
ln -s ${sem_source_dir}/DATA/* ${sem_build_dir}/DATA/

# copy setup directory
cp -r ${sem_source_dir}/setup ${sem_build_dir}/
chmod u+w -R setup

#====== use sem_config_dir
cd $sem_build_dir/DATA
rm -f Par_file CMTSOLUTION STATIONS
cp $sem_data_dir/Par_file .
cp $sem_data_dir/CMTSOLUTION .
cp $sem_data_dir/FORCESOLUTION .
cp $sem_data_dir/STATIONS .

cd $sem_build_dir/setup
rm -f constants.h.in
cp $sem_setup_dir/constants.h.in .

#====== build SEM
cd $sem_build_dir

#--stampede2
#NOTE: use -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 will cause 'double free' error. The code works with '-xCORE-AVX2'.
#./configure FC=ifort MPIFC=mpiifort CC=icc MPI_INC="${TACC_IMPI_INC}" FCFLAGS="-xCORE-AVX2" CFLAGS="-xCORE-AVX2"
#./configure FC=ifort MPIFC=mpiifort CC=icc MPI_INC="${TACC_IMPI_INC}" FCFLAGS="-xCORE-AVX2 -traceback" CFLAGS="-xCORE-AVX2 -traceback"
#./configure FC=mpif90 MPIFC=mpif90 CC=mpicc # FCFLAGS="-xCORE-AVX2 -O3" CFLAGS="-xCORE-AVX2 -O3"
# HIP_PATH="/public/software/compiler/dtk-24.04"
# HDF5_PATH=${HDF5}
# HDF5_PATH=/public/software/mathlib/hdf5/1.12.2/gnu_fortran_parallel
# ./configure --with-hip=MI50 --with-hdf5 FC=mpif90 MPIFC=mpif90 CC=mpicc FCFLAGS="-O3" CFLAGS="-O3" HIP_FLAGS="-fPIC -O2" HDF5_INC="${HDF5_PATH}/include" HDF5_LIBS="-L${HDF5_PATH}/lib" # -L${HIP_PATH}/lib"
# ./configure --with-hip=MI50 FC=mpif90 MPIFC=mpif90 CC=mpicc FCFLAGS="-O3" CFLAGS="-O3" HIP_FLAGS="-fPIC -O2"
if [ "$processor_type" == cpu ]
then
  echo "[INFO] Compile with CPU"
  ./configure FC=mpif90 MPIFC=mpif90 CC=mpicc FCFLAGS="-O3" CFLAGS="-O3"
elif [ "$processor_type" == dcu ]
then
  echo "[INFO] Compile with DCU"
  ./configure --with-hip=MI50 FC=mpif90 MPIFC=mpif90 CC=mpicc FCFLAGS="-O3" CFLAGS="-O3" HIP_FLAGS="-fPIC -O2"
else
  echo "[ERROR] Unknown processor type: "$processor_type""
  exit -1
fi

# make clean
# make all  > >(tee make.out) 2> >(tee make.err >&2)
make xmeshfem3D xspecfem3D xcombine_vol_data_vtu > >(tee make.out) 2> >(tee make.err >&2)

#END

