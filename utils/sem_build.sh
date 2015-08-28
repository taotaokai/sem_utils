#!/bin/bash

wkdir=$(pwd)

data_dir=DATA
build_dir=specfem3d_globe

# prepare data_dir
if [ ! -f $data_dir/Par_file ] || \
   [ ! -f $data_dir/CMTSOLUTION ] || \
   [ ! -f $data_dir/STATIONS ]
then
  echo "[ERROR] missing one or more control files in $data_dir"
  exit -1
fi

if [ ! -d DATABASES_MPI ] || [ ! -d OUTPUT_FILES ]
then
  mkdir DATABASES_MPI OUTPUT_FILES
fi

# prepare build_dir
cd $build_dir/

rm -rf DATABASES_MPI OUTPUT_FILES
ln -s $wkdir/DATABASES_MPI .
ln -s $wkdir/OUTPUT_FILES .

cd DATA
ln -sf $wkdir/DATA/Par_file .
ln -sf $wkdir/DATA/CMTSOLUTION .
ln -sf $wkdir/DATA/STATIONS .

# build
cd $wkdir/$build_dir

./configure FC=mpif90 MPIFC=mpif90 FCFLAGS="-O3"
make clean
make xcreate_header_file xmeshfem3D xspecfem3D xcombine_vol_data_vtk |\
  tee make.log