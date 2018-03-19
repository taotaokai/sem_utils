#!/bin/bash

wkdir=$(pwd)

control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}
dep_jobs=${3:--1}

# source control_file
source $control_file

#====== make slurm jobs
#work_flow=green,misfit,srcfrechet,dgreen,search
#work_flow=green
#work_flow=misfit
#work_flow=srcfrechet
#work_flow=dgreen
#work_flow=search

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  event_dir=$wkdir/$event_id

  chmod u+w -R $event_dir/DATA
  mkdir -p $event_dir/DATA

  # copy initial CMTSOLUTION file
  cmt_file=$source_dir/${event_id}.cmt
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]; then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi
  cp $cmt_file $event_dir/DATA/CMTSOLUTION.init

  # copy misfit_par file
  cp $misfit_par_dir/${event_id}_misfit_par.py $event_dir/DATA/misfit_par.py

  # create batch scripts
  #$utils_dir/make_source_iteration.sh $event_id
  #$wkdir/make_source_iteration.sh $event_id
  $sem_utils_dir/source_inversion/make_slurm_jobs_for_source_inversion.sh $control_file $event_id

done