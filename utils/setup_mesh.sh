#!/bin/bash

# setup mesh folders, generate the batch script to run SEM meshfem3D

#====== command line args
control_file=${1:?[arg] need control_file}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: " $control_file
    exit 1
fi
control_file=$(readlink -f $control_file)

# load parameters from control_file
source ${control_file}

#====== create mesh dir
if [ ! -d "${mesh_dir}" ]
then
    mkdir -p $mesh_dir
else
    echo "[ERROR] mesh_dir=$mesh_dir already exits!"
    exit 1
fi

#====== setup mesh dir
cd $mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

# create necessary files in mesh/DATA/
cd $mesh_dir/DATA
# link data files: topography, bathymetry, etc.
ln -sf $sem_build_dir/DATA/* ./

# modify Par_file
rm Par_file
cp -L $sem_config_dir/DATA/Par_file .
# SAVE_MESH_FILES = .false.
sed -i "/^[\s]*SAVE_MESH_FILES/s/=.*/= .false./" Par_file
# MODEL = GLL
sed -i "/^[\s]*MODEL/s/=.*/= GLL/" Par_file

# link previous model directory to DATA/GLL
rm $mesh_dir/DATA/GLL
ln -sf $prev_model_dir $mesh_dir/DATA/GLL

# backup Par_file into OUTPUT_FILES/
cp -L $mesh_dir/DATA/Par_file $mesh_dir/OUTPUT_FILES

#END
