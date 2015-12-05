#!/bin/bash

# create build directory

source_dir=${1:?[arg] need source folder}
#source_dir=~/tools/specfem3d_globe_modified
build_dir=${2:-specfem3d_globe}

source_dir=$(readlink -f $source_dir)
build_dir=$(readlink -f $build_dir)

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