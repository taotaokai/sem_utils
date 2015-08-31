#!/bin/bash

# create build directory

wkdir=$(pwd)

#source_dir=~/research/waveform_tomography/codes/specfem3d_globe_modified
source_dir=~/tools/specfem3d_globe_modified

build_dir=$wkdir/${1:-specfem3d_globe}

# prepare build_dir
rm -rf $build_dir
mkdir $build_dir
cd $build_dir

ln -s ${source_dir}/src ./
ln -s ${source_dir}/config* ./
ln -s ${source_dir}/flags.guess ./
ln -s ${source_dir}/install-sh ./
ln -s ${source_dir}/Makefile.in ./

mkdir DATA; cd DATA
ln -s ${source_dir}/DATA/* ./
cd ../

# copy setup directory
cp -r ${source_dir}/setup ./