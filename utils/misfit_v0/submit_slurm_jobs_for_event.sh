#!/bin/bash

# submit sequential slurm jobs for each event

wkdir=$(pwd)
utils_dir=~/seiscode/sem_utils/utils/misfit_v0

event_list=${1:?[arg]need event_list}
job_names=${2:?[arg]need job names, comma seperated (e.g. syn,misfit,kernel)}
job_dep=${3:--1} # dependent job ID's

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"
  event_dir=$wkdir/$event_id
  slurm_dir=$event_dir/slurm

  $utils_dir/submit_slurm_jobs.sh $slurm_dir $job_names $job_dep 

done
