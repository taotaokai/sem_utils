#!/bin/bash

# plot xsections for model relative difference
sem_utils=~/seiscode/sem_utils

model_dir=${1:?[arg]need model_dir (for *.nc)}
slice_list=${2:?[arg]need slice_list}
title=${3:?[arg]need title}
out_dir=${4:?[arg]need out_dir}

awk 'NF&&$1!~/#/' $slice_list |\
while read lat0 lon0 azimuth theta0 theta1 ntheta r0 r1 nr nc_tag 
do
  echo
  echo "#====== $nc_tag: $lat0 $lon0 $azimuth $theta0 $theta1 $ntheta $r0 $r1 $nr"
  echo
  $sem_utils/utils/xsection/plot_slice_gcircle_dlnvs_kappa_thomsen_elliptic.py \
    $model_dir $nc_tag \
    $lat0 $lon0 $azimuth \
    "$title" $out_dir
done