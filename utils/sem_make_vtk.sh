#!/bin/bash

#Usage: xcombine_vol_data slice_list filename input_topo_dir input_file_dir 
#         output_dir high/low-resolution [region]

if [ $# -eq 0 ]
then
cat<<EOF
Usage: $0 <model_tags> <nproc> <mesh_dir> <model_dir> <out_dir>
EOF
  exit
fi

model_name_array=${1:-vp,vs}
nproc=${2:-1}
mesh_dir=${3:-DATABASES_MPI}
model_dir=${4:-DATABASES_MPI}
vtk_dir=${5:-VTK}

sem_bin=specfem3d_globe/bin

if [ ! -d "$vtk_dir" ]; then
  mkdir $vtk_dir 
fi

#nproc=$(grep NPROCTOT_VAL OUTPUT_FILES/values_from_mesher.h | awk '{print $NF}')
seq 0 $(($nproc - 1)) > $vtk_dir/slice_list

for model_name in $(echo $model_name_array | sed "s/,/ /g")
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