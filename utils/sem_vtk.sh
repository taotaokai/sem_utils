#!/bin/bash

# make vtk file of SEM GLL model

#====== command line args
model_names=${1:?[arg] need model_name (e.g. mu_kernel,lamda_kernel,rho_kernel)}
sem_dir=${2:?[arg] need sem_dir (for bin/xcombine_vol_data_vtk)}

topo_dir=DATABASES_MPI
model_dir=DATABASES_MPI
out_dir=vtk
iresolution=0
iregion=1

if [ ! -d "$out_dir" ]
then
  mkdir $out_dir
fi

# make slice_list
nproc=$(ls $topo_dir/*_reg1_solver_data.bin | wc -l)
seq 0 $((nproc - 1)) > $out_dir/slice.list

#for tag in betav_kernel #betah_kernel alphav_kernel alphah_kernel
for tag in ${model_names//,/ }
do
    $sem_dir/bin/xcombine_vol_data_vtk \
        $out_dir/slice.list $tag $topo_dir $model_dir \
        $out_dir $iresolution $iregion
done
