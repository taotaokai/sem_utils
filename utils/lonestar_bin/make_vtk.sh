#!/bin/bash

# make vtk 
# for forward simulation

#====== command line args
control_file=${1:?must provide control_file}
evid=${2:?must provide evid}
model_names=${3:-mu_kernel,lamda_kernel,rho_kernel}

# source parameters in control_file
source ${control_file}

topo_dir=${mesh_dir}/DATABASES_MPI
model_dir=DATABASES_MPI
out_dir=vtk
iresolution=0
iregion=1

#====== 
event_dir=${iter_dir}/$evid

cd $event_dir
if [ ! -d $out_dir ]
then
    mkdir $out_dir 
fi

# make slice_list
seq 0 $((nproc - 1)) > $out_dir/slice.list
#awk '$1=="NPROC_XI"{a=$3};
#    $1=="NPROC_ETA"{b=$3};
#    END{n=a*b; for(i=0;i<n;i++) print i}' DATA/Par_file \

#for tag in betav_kernel #betah_kernel alphav_kernel alphah_kernel
for tag in ${model_names//,/ }
do
    $build_dir/bin/xcombine_vol_data_vtk \
        $out_dir/slice.list $tag $topo_dir $model_dir \
         $out_dir $iresolution $iregion
done
