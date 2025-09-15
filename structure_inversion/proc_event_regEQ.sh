#!/bin/bash

# structure inversion  for regional earthquake

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
control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}

# load parameters in control_file
source $control_file

if [ ! -d "$updated_model_dir" ]
then
  mkdir -p $updated_model_dir
fi

#------ process each event
# for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  echo "====== $event_id"

  # create event dir
  event_dir=$iter_dir/events/$event_id
  sem_data_dir=$event_dir/DATA
  if [ -d "$sem_data_dir" ]
  then
    chmod u+w -R $event_dir/DATA
  fi
  mkdir -p $event_dir/DATA

  mkdir -p $event_dir

  # copy CMTSOLUTION file
  cmt_file=$(find -L $source_dir -path "*/iter??/events/${event_id}/misfit/CMTSOLUTION.updated" | sort | tail -n1)
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi

  chmod u+w $event_dir/DATA/CMTSOLUTION.init
  cp $cmt_file $event_dir/DATA/CMTSOLUTION.init
  # chmod a-w $event_dir/DATA/CMTSOLUTION.init

  # copy STATIONS
  station_file=$data_dir/$event_id/STATIONS
  if [ ! -f "$station_file" ]; then
    echo "[ERROR] $station_file not found"
    exit -1
  fi
  cp $station_file $event_dir/DATA/STATIONS

  # copy misfit_par file
  # cp $misfit_par_dir/${event_id}_misfit.yaml $event_dir/DATA/misfit.yaml
  cp $misfit_par_dir/misfit.yaml $event_dir/DATA/misfit.yaml

  # create batch scripts
  $sem_utils_dir/structure_inversion/make_slurm_jobs_for_event_regEQ.sh $control_file $event_id

done
