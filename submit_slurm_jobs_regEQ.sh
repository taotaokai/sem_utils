#!/bin/bash

# submit sequential slurm jobs for each event

wkdir=$(pwd)

# utils_dir=$(readlink -f ${0} | sed 's/\/[^\/]*$/\//')

control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}
job_names=${3:?[arg]need job names, comma seperated (e.g. green,misfit,srcfrechet,dgreen,search)}
job_dep=${4:--1} # dependent job ID's

# load parameters in control_file
source $control_file

# for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
for event_id in $(awk 'NF&&$1!~/#/' $event_list)
do

  echo "====== $event_id"
  event_dir=$SEM_iter_dir/events/$event_id
  slurm_dir=$event_dir/slurm

  $SEM_utils_dir/submit_slurm_jobs_sequential.sh $slurm_dir $job_names $job_dep

done
