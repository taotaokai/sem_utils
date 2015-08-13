#!/bin/bash

#Usage: xcombine_vol_data slice_list filename input_topo_dir input_file_dir 
#         output_dir high/low-resolution [region]

model_name_array=${1:-vp,vs}
syn_dir=${2:-.}
model_dir=$syn_dir/${3:-DATABASES_MPI}
mesh_dir=$syn_dir/${4:-DATABASES_MPI}
vtk_dir=$syn_dir/${5:-VTK}

sem_bin=specfem3d_globe/bin

echo "### syn_dir=$syn_dir vtk_dir=$vtk_dir model_name=$model_name_array"

if [ ! -d "$vtk_dir" ]; then
  mkdir $vtk_dir 
fi

nproc=$(grep NPROCTOT_VAL OUTPUT_FILES/values_from_mesher.h | awk '{print $NF}')
seq 0 $(($nproc - 1)) > $vtk_dir/slice_list

for model_name in $(echo $model_name_array | sed "s/,/ /")
do
  echo "#========================== "
  echo "# process $model_name"
  echo "#========================== "
  $sem_bin/xcombine_vol_data_vtk \
    $vtk_dir/slice_list \
    $model_name \
    $mesh_dir \
    $model_dir \
    $vtk_dir 0 1
done

#bin/xcombine_vol_data_vtk \
#  VTK/slice_list vs DATABASES_MPI MODEL_SLAB VTK 1 1