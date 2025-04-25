#!/bin/bash

# source inversion for regional earthquake

# Work flows contain:
#work_flow=green,misfit,srcfrechet,dgreen,search
#work_flow=green
#work_flow=misfit
#work_flow=srcfrechet
#work_flow=dgreen
#work_flow=search

#------ read command line args
control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}

# load parameters in control_file
source $control_file

#------ process each event
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  echo "====== $event_id"

  # create event dir
  event_dir=$iter_dir/$event_id
  chmod u+w -R $event_dir/DATA
  mkdir -p $event_dir/DATA

  # copy initial CMTSOLUTION file
  cmt_file=$source_dir/${event_id}.cmt
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]; then
    echo "[WARN] $cmt_file does NOT exist!"
    exit -1
  fi
  rm $event_dir/DATA/CMTSOLUTION.init
  cp -L $cmt_file $event_dir/DATA/CMTSOLUTION.init

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
  $sem_utils_dir/source_inversion_regEQ/make_slurm_jobs_for_event_regEQ.sh $control_file $event_id

done
