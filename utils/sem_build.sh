#!/bin/bash

data_dir=${1:-DATA}
build_dir=${2:-specfem3d_globe}

echo "### data_dir=$data_dir build_dir=$build_dir"

data_dir=$(readlink -f $data_dir)
build_dir=$(readlink -f $build_dir)

# prepare data_dir
if [ ! -f $data_dir/Par_file ] || \
   [ ! -f $data_dir/CMTSOLUTION ] || \
   [ ! -f $data_dir/STATIONS ]; then
  echo "[ERROR] missing one or more control files in $data_dir/"
  exit -1
fi

# prepare build_dir
cd $build_dir

rm -rf DATA OUTPUT_FILES
mkdir DATA OUTPUT_FILES

cd DATA
ln -s $data_dir/* .

# build
cd $build_dir

./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-O3"
make clean
make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk |\
  tee make.log