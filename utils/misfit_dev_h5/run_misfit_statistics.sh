#!/bin/bash

work_dir=${1:?[arg] work_dir, e.g. ./}
inv_type=${2:?[arg] inv_type, e.g. source}
stage_num=${3:?[arg] stage_num, e.g. 0}
max_niter=${4:?[arg] max_niter, e.g. 0}
event_list=${5:?[arg] events.lst, e.g. events.lst}

num=$(echo $stage_num | awk '{printf "%02d",$1}')
stage_dir=${work_dir}/stage${num}_${inv_type}
for iter_num in $(seq 0 $max_niter)
do
  num=$(echo $iter_num | awk '{printf "%02d",$1}')
  iter_dir=$stage_dir/iter${num}

  for event_name in $(awk '$1!~/^#/{print $1}' $event_list)
  do
    misfit_dir=$iter_dir/events/$event_name/misfit
    if [ ! -f $misfit_dir/misfit.h5 ]; then
      continue
    fi
    echo python sem_utils/misfit/make_misfit_statistics.py \
      $misfit_dir/misfit.h5 \
      $misfit_dir/misfit.csv \
      --stage $stage_num --iter $iter_num --type $inv_type
  done
done
