#!/bin/bash

#Usage: xcombine_vol_data slice_list filename input_topo_dir input_file_dir 
#         output_dir high/low-resolution [region]

syn_dir=${1:-.}
model_name_array=${2:-vp,vs}
sem_bin=specfem3d_globe/bin

vtk_dir=$syn_dir/vtk

echo "### syn_dir=$syn_dir vtk_dir=$vtk_dir model_name=$model_name_array"

if [ ! -d "$vtk_dir" ]; then
  mkdir $syn_dir/VTK
fi

nproc=$(grep NPROCTOT_VAL $syn_dir/OUTPUT_FILES/values_from_mesher.h | awk '{print $NF}')
seq 0 $(($nproc - 1)) > $vtk_dir/slice_list

for model_name in $(echo $model_name_array | sed "s/,/ /")
do
  echo "#========================== "
  echo "# process $model_name"
  echo "#========================== "
  $sem_bin/xcombine_vol_data_vtk \
    $vtk_dir/slice_list \
    $model_name \
    $syn_dir/DATABASES_MPI \
    $syn_dir/DATABASES_MPI \
    $vtk_dir 0 1
done

#bin/xcombine_vol_data_vtk \
#  VTK/slice_list vs DATABASES_MPI MODEL_SLAB VTK 1 1