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

if [ ! -d "$updated_cmt_dir" ]
then
  mkdir -p $updated_cmt_dir
fi

#------ process each event
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

  # copy Par_file
  par_file=${sem_config_dir}/DATA/Par_file
  if [ ! -f "$par_file" ]
  then
    echo "[ERROR] $par_file does NOT exist!"
    exit -1
  fi
  cp $par_file $event_dir/DATA/Par_file

  # copy initial CMTSOLUTION file
  if [ x${iter_num} == x ]
  then
    echo "[ERROR] iter_num not set!"
    exit -1
  fi
  if [ "$iter_num" -eq 0 ]
  then
    cmt_file=$data_dir/${event_id}/CMTSOLUTION.ecef
  else
    cmt_file=$initial_cmt_dir/${event_id}.cmt
  fi
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file does NOT exist!"
    exit -1
  fi

  init_cmt_file=$event_dir/DATA/CMTSOLUTION.init
  if [ -f "$init_cmt_file" ]
  then
    rm -f $init_cmt_file
  fi
  cp -L $cmt_file $init_cmt_file

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
