#!/bin/bash

dmodel=${1:?[arg]need model types (e.g. perturb,random)}

wkdir=$(pwd)

sem_utils=~/seiscode/sem_utils
#mpiexec=ibrun
mpiexec="mpirun -np 256"

base_dir=~/NAmer
par_dir=$base_dir/sem_config/DATA
sem_dir=$base_dir/specfem3d_globe

#for dmodel in dvp dvsv dvsh
#for dmodel in perturb
for dmodel in ${dmodel//,/ }
do
  echo ====== $dmodel

  mesh_dir=$wkdir/mesh_${dmodel}
  model_dir=$wkdir/model_${dmodel}

  if [ ! -d "$model_dir" ]
  then
    echo "[ERROR] $model_dir does not exist!"
    exit -1
  fi

  rm -rf $mesh_dir

  $sem_utils/utils/sem_mesh.sh $mesh_dir $par_dir $sem_dir "$mpiexec"

  sed -i "/^MODEL/s/=.*/= GLL/" $mesh_dir/DATA/Par_file
  ln -s $model_dir $mesh_dir/DATA/GLL

done