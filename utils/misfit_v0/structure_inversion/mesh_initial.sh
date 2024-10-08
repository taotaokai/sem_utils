#!/bin/bash

base_dir=~/NAmer

wkdir=$(pwd)

sem_utils=~/seiscode/sem_utils
par_dir=$base_dir/sem_config/DATA
sem_dir=$base_dir/specfem3d_globe

echo ====== initial

mesh_dir=$wkdir/mesh
model_dir=$wkdir/model

if [ ! -d "$model_dir" ]
then
  echo "[ERROR] $model_dir does not exist!"
  exit -1
fi

rm -rf $mesh_dir

$sem_utils/utils/sem_mesh.sh $mesh_dir $par_dir $sem_dir "mpirun -np 256"

ln -s $model_dir $mesh_dir/DATA/GLL