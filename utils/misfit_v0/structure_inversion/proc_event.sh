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
for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do
  echo "====== $event_id"

  event_dir=$iter_dir/$event_id
  mkdir -p $event_dir

  #rm -rf $event_dir
  chmod u+w -R $event_dir/DATA
  mkdir -p $event_dir/DATA

  # copy CMTSOLUTION file
  #cmt_file=$(find -L $source_dir -path "*/iter0?/CMTSOLUTION_updated/${event_id}.cmt" | sort | tail -n1)
  cmt_file=$(ls $source_dir/iter0?/CMTSOLUTION_updated/${event_id}.cmt | sort | tail -n1)
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]; then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi
  chmod u+w $event_dir/DATA/CMTSOLUTION.init
  cp $cmt_file $event_dir/DATA/CMTSOLUTION.init

  # copy STATIONS
  station_file=$data_dir/$event_id/data/STATIONS
  if [ ! -f "$station_file" ]; then
    echo "[ERROR] $station_file not found"
    exit -1
  fi
  cp $station_file $event_dir/DATA/STATIONS
  #remove stations too close to the mesh western boundary
  #awk '$1!~/#/{if($3<40 && $4>90) print $0; if($3>=40 && $4>=87) print $0;}' $station_file > $event_dir/DATA/STATIONS

  # copy misfit_par file
  cp $misfit_par_dir/${event_id}_misfit_par.py $event_dir/DATA/misfit_par.py

  # create batch scripts
  $sem_utils_dir/utils/misfit_v0/structure_inversion/make_slurm_jobs_for_event.sh $control_file $event_id

done