#!/bin/bash

sem_utils=~/seiscode/sem_utils

nproc=256
#nproc=1

mesh_dir=mesh_spherical
model_dir=model
#model_dir=model_new_starting

out_dir=model_radial_profiles
mkdir $out_dir

for iproc in $(seq -f "%03.0f" 0 $((nproc - 1)) )
do

  echo ====== iproc = $iproc
  #$sem_utils/bin/xsem_get_radius $iproc $mesh_dir $model_dir vsv,vpv,eps,gamma > $out_dir/model_proc000${iproc}.txt
  $sem_utils/bin/xsem_get_radius $iproc $mesh_dir $model_dir vsv,vpv,dlnvs,kappa,eps,gamma > $out_dir/model_proc000${iproc}.txt

done