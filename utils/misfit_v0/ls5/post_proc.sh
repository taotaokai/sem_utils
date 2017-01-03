#!/bin/bash

# submit kernel post-process jobs files 

wkdir=$(pwd)

job_id=${1:--1}

#for job in dmodel_threshold dmodel_smooth dmodel_scale dmodel_statis
#for job in kernel_threshold kernel_smooth pcg_dmodel dmodel_scale kernel_statis model_perturb
for job in kernel_sum kernel_threshold kernel_smooth pcg_dmodel dmodel_scale kernel_statis model_perturb
do 

  echo
  echo "====== job: $job"
  echo

  job_file=$job.job

  if [ ! -f $job_file ]
  then
    echo "[ERROR] job file ($job_file) does NOT exist!"
    exit -1
  fi

  log_file=$job_file.submit

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
