#!/bin/bash

# submit sequential slurm jobs

slurm_dir=${1:?[arg]need slurm_dir (contain *.job)}
job_names=${2:?[arg]need job names, comma seperated (e.g. syn,misfit,kernel)}
job_dep=${3:--1} # dependent job ID's

job_id=$job_dep
for work in ${job_names//,/ }
do 

  echo
  echo "====== job: $work"

  job_file=$slurm_dir/${work}.job

  if [ ! -f $job_file ];then
    echo "[ERROR] job file ($job_file) does NOT exist!"
    exit -1
  fi

  log_file=${job_file}.submit

  sbatch --dependency=afterok:$job_id $job_file | tee $log_file

  grep -i error $log_file > /dev/null
  if [ $? -eq 0 ]
  then
    echo "[ERROR] failed to submit $job_file"
    exit -1
  fi

  job_id=$(grep "Submitted batch job" $log_file | awk '{print $NF}')
  echo "job_id=$job_id"

done
