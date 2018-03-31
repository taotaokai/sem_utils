#!/bin/bash

sem_utils=~/seiscode/sem_utils

nproc=144
#nproc=1

mesh_dir=mesh_EARA2014_spherical/DATABASES_MPI
model_dir=model_EARA2014

out_dir=radial_model
mkdir $out_dir

for iproc in $(seq -f "%03.0f" 0 $((nproc - 1)) )
do

  echo ====== iproc = $iproc
  #$sem_utils/bin/xsem_get_radius $iproc $mesh_dir $model_dir vsv,vpv,eps,gamma > $out_dir/model_proc000${iproc}.txt
  $sem_utils/bin/xsem_get_radius $iproc $mesh_dir $model_dir vsv,vsh,vpv,vph > $out_dir/model_proc000${iproc}.txt

done