#!/bin/bash

# submit sequential slurm jobs for each event

wkdir=$(pwd)

utils_dir=$(readlink -f ${0} | sed 's/\/[^\/]*$/\//') #~/seiscode/specfem3d_globe/utils_ktao

event_list=${1:?[arg]need event_list}
job_names=${2:?[arg]need job names, comma seperated (e.g. syn,misfit,kernel)}
job_dep=${3:--1} # dependent job ID's

for event_id in $(awk -F"|" 'NF&&$1!~/#/{printf "%s.%s.%s.%s\n", $1,$2,$3,$4}' $event_list)
do

  echo "====== $event_id"
  event_dir=$wkdir/$event_id
  slurm_dir=$event_dir/slurm

  $utils_dir/submit_slurm_jobs_sequential.sh $slurm_dir $job_names $job_dep 

done
