#!/bin/bash

# plot xsections for model relative difference
sem_utils=~/seiscode/sem_utils

model_dir=${1:?[arg]need model_dir (*.nc)}
slice_list=${2:?[arg]need slice_list}
title=${3:?[arg]need title}
out_dir=${4:?[arg]need out_dir}

awk 'NF&&$1!~/#/' $slice_list |\
while read lat0 lat1 nlat lon0 lon1 nlon depth nc_tag 
do

  echo
  echo "#====== $nc_tag: $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth"
  echo

  $sem_utils/utils/xsection/plot_slice_sphere.py \
    $model_dir \
    $nc_tag \
    "$title" $out_dir
done