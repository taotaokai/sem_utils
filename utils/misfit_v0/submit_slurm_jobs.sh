#!/bin/bash
# submit squential jobs files for slurm 

wkdir=$(pwd)

event_id=${1:?[arg]need event_id}
work_flow=${2:?[arg]need work_flow (e.g. syn,misfit,kernel,hess)}
job_dep=${3:--1} # dependent jobs

event_dir=$wkdir/$event_id

job_id=$job_dep
for work in ${work_flow//,/ }
do 

  echo
  echo "====== job: $work"

  job_file=$event_dir/slurm/$work.job

  if [ ! -f $job_file ];then
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
