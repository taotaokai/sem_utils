#!/bin/bash

# structure inversion 

# Work flows contain:
#--  kernel
#work_flow=syn,misfit,kernel
#work_flow=syn
#work_flow=misfit
#work_flow=kernel
#-- line search
#work_flow=perturb
#work_flow=search
#-- hessian-random model product
#work_flow=perturb_random,misfit_random,kernel_random
#work_flow=perturb_random
#work_flow=misfit_random
#work_flow=kernel_random

#------ read command line args
control_file=${1:?[arg] need control_file}
event_list=${2:?[arg]need event_list}

#------ source control_file
source $control_file

#------ process each event
#for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
for event_id in $(awk -F"|" 'NF&&$1!~/#/{printf "%s.%s.%s.%s\n", $1,$2,$3,$4}' $event_list)
do
  echo "====== $event_id"

  event_dir=$iter_dir/$event_id
  mkdir -p $event_dir

  #rm -rf $event_dir
  chmod u+w -R $event_dir/DATA
  mkdir -p $event_dir/DATA

  # copy CMT and FORCESOLUTION file
  cmt_file=$(ls ${data_dir}/${event_id}/data/CMTSOLUTION)
  force_file=$(ls ${data_dir}/${event_id}/data/FORCESOLUTION)
  echo ------ use: $(readlink -f $force_file $cmt_file)
  if [ ! -f "$cmt_file" ] || [ ! -f "$force_file" ] ; then
    echo "[ERROR] $cmt_file and/or $force_file not found"
    exit -1
  fi
  cp $cmt_file $event_dir/DATA/CMTSOLUTION
  cp $force_file $event_dir/DATA/FORCESOLUTION

  # copy STATIONS
  station_file=$data_dir/$event_id/data/STATIONS
  if [ ! -f "$station_file" ]; then
    echo "[ERROR] $station_file not found"
    exit -1
  fi
  cp $station_file $event_dir/DATA/STATIONS

  # copy misfit_par file
  cp $misfit_par_dir/${event_id}_misfit_par.py $event_dir/DATA/misfit_par.py

  # create batch scripts
  $sem_utils_dir/structure_inversion/make_slurm_jobs_for_event_noiseCC.sh $control_file $event_id

done