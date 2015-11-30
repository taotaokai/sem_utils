#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:-control_file}
evid=${2:?must provide event id}

# source parameters in control_file
source ${control_file}

echo "Start model updating [$(date)]."

#====== loop each event
event_dir=${iter_dir}/$evid

#-- create kernel mask file
cd $event_dir
$sem_utils/bin/xsem_make_kernel_mask \
    $nproc DATABASES_MPI "OUTPUT_forward/source.vtk" DATABASES_MPI \
    $source_mask_radius $stop_depth $pass_depth

#-- reduce cijkl kernel to lamda,mu kernel
cd $event_dir
$sem_utils/bin/xsem_kernel_cijkl_to_lamda_mu \
    $nproc DATABASES_MPI DATABASES_MPI DATABASES_MPI

#-- sum event kernels with kernel masks
cd ${iter_dir}
#mkdir kernel_sum
# do something

#-- get model udpate direction
cd ${iter_dir}
mkdir model_update_direction
# specially for 1 event
ln -sf $evid/DATABASES_MPI kernel_sum
$sem_utils/bin/xsem_kernel_lamda_mu_to_dmodel_sd \
    $nproc $mesh_dir/DATABASES_MPI kernel_sum \
    model_update_direction 1.0  1

#-- make new model
cd ${iter_dir}
mkdir model_new
$sem_utils/bin/xsem_add_dmodel_lamda_mu_to_tiso \
    $nproc $mesh_dir/DATABASES_MPI $mesh_dir/DATA/GLL \
    model_update_direction model_new \
    $max_dlnv_allowed $force_max_dlnv_allowed

#
echo "The model update is finished [$(date)]."
