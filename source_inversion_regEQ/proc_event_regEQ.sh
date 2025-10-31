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
if [ ! -f "$control_file" ]
then
  echo "[ERROR] $control_file NOT found!"
  exit -1
fi
source $control_file
if [ $? -ne 0 ]
then
  echo "[ERROR] source $control_file failed!"
  exit -1
fi

if [ ! -d "$SEM_iter_dir/CMTSOLUTION_initial" ]
then
  mkdir -p $SEM_iter_dir/CMTSOLUTION_initial
fi
if [ ! -d "$SEM_iter_dir/CMTSOLUTION_updated" ]
then
  mkdir -p $SEM_iter_dir/CMTSOLUTION_updated
fi

#------ process each event
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  echo "====== $event_id"

  # create event dir
  event_dir=$SEM_iter_dir/events/$event_id
  if [ -d "$event_dir/DATA" ]
  then
    chmod u+w -R $event_dir/DATA # make DATA/ writable
  else
    mkdir -p $event_dir/DATA
  fi

  # copy Par_file
  par_file=${SEM_config_dir}/DATA/Par_file
  # par_file=${SEM_mesh_dir}/DATA/Par_file
  if [ ! -f "$par_file" ]
  then
    echo "[ERROR] $par_file NOT found!"
    exit -1
  fi
  cp -L $par_file $event_dir/DATA/Par_file

  # copy initial CMTSOLUTION file
  if [ x${SEM_iter_num} == x ]
  then
    echo "[ERROR] iter_num NOT set!"
    exit -1
  fi
  if [ "$SEM_iter_num" -eq 0 ]
  then
    cmt_file=$SEM_data_dir/${event_id}/CMTSOLUTION.ecef
  else
    cmt_file=$SEM_prev_iter_dir/CMTSOLUTION_updated/${event_id}.cmt
  fi
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file NOT found!"
    exit -1
  fi
  # init_cmt_file=$event_dir/DATA/CMTSOLUTION.init
  # if [ -f "$init_cmt_file" ]
  # then
  #   rm -f $init_cmt_file
  # fi
  cp -L $cmt_file $SEM_iter_dir/CMTSOLUTION_initial/${event_id}.cmt

  echo ------ use: $(readlink -f $cmt_file)

  # copy STATIONS
  station_file=$SEM_data_dir/$event_id/STATIONS
  if [ ! -f "$station_file" ]; then
    echo "[ERROR] $station_file NOT found!"
    exit -1
  fi
  # if [ -f "$event_dir/DATA/STATIONS" ]
  # then
  #   chmod 644 $event_dir/DATA/STATIONS
  # fi
  cp -L $station_file $event_dir/DATA/STATIONS

  # copy misfit_par file
  # cp $misfit_par_dir/${event_id}_misfit.yaml $event_dir/DATA/misfit.yaml
  if [ ! -f "$SEM_misfit_par_dir/misfit.yaml" ]; then
    echo "[ERROR] $SEM_misfit_par_dir/misfit.yaml NOT found!"
    exit -1
  fi
  cp -L $SEM_misfit_par_dir/misfit.yaml $event_dir/DATA/misfit.yaml

  # create batch scripts
  $SEM_utils_dir/source_inversion_regEQ/make_slurm_jobs_for_event_regEQ.sh $control_file $event_id

done