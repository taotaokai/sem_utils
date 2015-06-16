#!/bin/bash

syn_dir=${1:-.}
build_dir=specfem3d_globe

echo "### syn_dir=$syn_dir build_dir=$build_dir"

syn_dir=$(readlink -f $syn_dir)
build_dir=$(readlink -f $build_dir)

# prepare syn_dir
if [ ! -f $syn_dir/DATA/Par_file ] || \
   [ ! -f $syn_dir/DATA/CMTSOLUTION ] || \
   [ ! -f $syn_dir/DATA/STATIONS ]; then
  echo "[ERROR] missing one or more control files in $syn_dir/DATA"
  exit -1
fi

if [ ! -d "$syn_dir/DATABASES_MPI" ]; then
  mkdir $syn_dir/DATABASES_MPI
fi
if [ ! -d "$syn_dir/OUTPUT_FILES" ]; then
  mkdir $syn_dir/OUTPUT_FILES
fi

# prepare build_dir
cd $build_dir

rm -rf DATABASES_MPI OUTPUT_FILES
ln -s $syn_dir/DATABASES_MPI ./
ln -s $syn_dir/OUTPUT_FILES ./

cd DATA
rm -f Par_file CMTSOLUTION STATIONS
ln -s $syn_dir/DATA/Par_file ./
ln -s $syn_dir/DATA/CMTSOLUTION ./
ln -s $syn_dir/DATA/STATIONS ./

# build 
cd $build_dir

make clean

./configure FC=mpif90 MPIFC=mpif90

make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk