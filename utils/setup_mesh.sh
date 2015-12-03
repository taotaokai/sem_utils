#!/bin/bash

# setup mesh folders, generate the batch script to run SEM meshfem3D

#====== command line args
control_file=${1:?must provide control_file}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: " $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)

# load parameters from control_file
source ${control_file}

#====== create mesh dir
if [ ! -d ${mesh_dir} ];then
    mkdir -p $mesh_dir
else
    echo "[WARNING] mesh_dir=$mesh_dir already exits!"
    exit -1
fi

#====== setup mesh dir
cd $mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

# create necessary files in mesh/DATA/
cd $mesh_dir/DATA
# all data files: topography, bathymetry, etc.
ln -sf $sem_build_dir/DATA/* ./
# modify Par_file
rm Par_file
cp -L $sem_config_dir/DATA/Par_file .
if [ ${iter} -eq 0 ] # save mesh file for starting model
then
    sed -i "/^SAVE_MESH_FILES/s/=.*/= .true./" Par_file
else
    sed -i "/^SAVE_MESH_FILES/s/=.*/= .false./" Par_file
fi

# link model directory to DATA/GLL
cd $mesh_dir/DATA
rm GLL
ln -sf $model_dir/DATABASES_MPI GLL

# backup parameter files into OUTPUT_FILES
cd $mesh_dir
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES

#END
