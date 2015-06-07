#!/bin/bash

#Usage: xcombine_vol_data slice_list filename input_topo_dir input_file_dir 
#         output_dir high/low-resolution [region]

mkdir VTK

nproc=$(grep NPROCTOT_VAL OUTPUT_FILES/values_from_mesher.h | awk '{print $NF}')
seq 0 $(($nproc - 1)) > VTK/slice_list

bin/xcombine_vol_data_vtk \
  VTK/slice_list vp DATABASES_MPI DATABASES_MPI VTK 0 1

#bin/xcombine_vol_data_vtk \
#  VTK/slice_list vs DATABASES_MPI MODEL_SLAB VTK 1 1