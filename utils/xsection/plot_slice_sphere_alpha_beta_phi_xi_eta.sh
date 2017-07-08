#!/bin/bash

# plot horizontal slices 
sem_utils=~/seiscode/sem_utils

nc_dir=${1:?[arg]need nc_dir (for *.nc)}
slice_list=${2:?[arg]need slice_list}
title=${3:?[arg]need title}
model_names=${4:?[arg]need model names (e.g. xi,beta,alpha)}
out_dir=${5:?[arg]need out_dir}

awk 'NF&&$1!~/#/' $slice_list |\
while read lat0 lat1 nlat lon0 lon1 nlon depth nc_tag 
do
  echo
  echo "#====== $nc_tag: $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth"
  echo

  out_fig=$out_dir/${nc_tag}.pdf

  $sem_utils/utils/xsection/plot_slice_sphere_alpha_beta_phi_xi_eta.py \
    $nc_dir/${nc_tag}.nc "$title ($depth km)" $out_fig
done