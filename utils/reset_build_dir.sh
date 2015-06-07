#!/bin/bash

code_dir=~/research/waveform_tomography/codes/specfem3d_globe_user_changes
build_dir=specfem3d_globe
wkdir=$(pwd)

# prepare code dir
cd $wkdir

rm -rf $build_dir bin
mkdir $build_dir
cd $build_dir

ln -s ${code_dir}/src ./
ln -s ${code_dir}/config* ./
ln -s ${code_dir}/flags.guess ./
ln -s ${code_dir}/install-sh ./
ln -s ${code_dir}/Makefile.in ./

cp -L -r ${code_dir}/setup ./ 
ln -s ../DATA ./
ln -s ../OUTPUT_FILES ./

# prepare DATA
cd $wkdir
mkdir DATA DATABASES_MPI OUTPUT_FILES
ln -s $build_dir/bin ./

cd DATA
ln -sf $code_dir/DATA/* ./