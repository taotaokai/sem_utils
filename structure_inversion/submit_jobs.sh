#!/bin/bash
set -x

control_file=${1:?[arg]need control_file, e.g. control_structure_dcu }
event_list=${2:?[arg]need event_list, e.g. events.lst }

source ${control_file}
if [ $? -ne 0 ]
then
  echo "[ERROR] failed to source $control_file"
  exit -1
fi

${SEM_utils_dir}/structure_inversion/proc_event_regEQ.sh ${control_file} ${event_list}

job_file=${SEM_iter_dir}/slurm/mesh.job
job_id=$(sbatch --parsable $job_file)
if [ $? -ne 0 ]
then
  echo "[ERROR] failed submitting $job_file"
  exit -1
fi
echo mesh ${job_id}  

job_list=submitted_jobs.lst
[ -f ${job_list} ] && rm ${job_list} || touch ${job_list}

dep_job=${job_id}
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  event_dir=$SEM_iter_dir/events/$event_id
  prev_job=${dep_job}
  for job_name in forward misfit kernel threshold
  do
    job_file=${event_dir}/slurm/${job_name}.job
    job_id=$(sbatch --parsable -d afterok:${prev_job} $job_file)
    if [ $? -ne 0 ]
    then
      echo "[ERROR] failed submitting $job_file"
      exit -1
    fi
    prev_job=${job_id}
    echo ${event_id} ${job_name} ${job_id}  
  done
  echo $job_id >> $job_list

  job_file=${event_dir}/slurm/plot.job
  job_id=$(sbatch --parsable -d afterok:${job_id} $job_file)
  if [ $? -ne 0 ]
  then
    echo "[ERROR] failed submitting $job_file"
    exit -1
  fi
  echo ${event_id} plot ${job_id}  

done

dep_jobs=$(awk '{printf ":%s", $1}' $job_list)
for job_name in kernel_sum kernel_precond dmodel model_perturb mesh_perturb
do
  job_file=${SEM_iter_dir}/slurm/${job_name}.job
  job_id=$(sbatch --parsable -d afterok${dep_jobs} $job_file)
  if [ $? -ne 0 ]
  then
    echo "[ERROR] failed submitting $job_file"
    exit -1
  fi
  dep_jobs=:$job_id
  echo ${job_name} ${job_id}  
done

job_list=submitted_job.lst
[ -f ${job_list} ] && rm ${job_list} || touch ${job_list}

for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  event_dir=$SEM_iter_dir/events/$event_id
  for job_name in perturb search
  do
    job_file=${event_dir}/slurm/${job_name}.job
    job_id=$(sbatch --parsable -d afterok${dep_jobs} $job_file)
    if [ $? -ne 0 ]
    then
      echo "[ERROR] failed submitting $job_file"
      exit -1
    fi
    dep_jobs=:$job_id
    echo ${event_id} ${job_name} ${job_id}  
  done
  echo $job_id >> $job_list
done

dep_jobs=$(awk '{printf ":%s", $1}' $job_list)
job_file=${SEM_iter_dir}/slurm/model_update.job
job_id=$(sbatch --parsable -d afterok${dep_jobs} $job_file)
if [ $? -ne 0 ]
then
  echo "[ERROR] failed to submit $job_file"
  exit -1
fi
echo model_update ${job_id}  

# ./sem_utils/submit_slurm_jobs_sequential.sh stage01_structure_dcu/iter01/slurm mesh
#
# ./sem_utils/submit_slurm_jobs_regEQ.sh control_file_structure_dcu events.lst forward,misfit,kernel,threshold,plot
#
# ./sem_utils/submit_slurm_jobs_sequential.sh stage01_structure_dcu/iter01/slurm kernel_sum,kernel_precond,dmodel,model_perturb,mesh_perturb
#
# ./sem_utils/submit_slurm_jobs_regEQ.sh control_file_structure_dcu events.lst perturb,search
#
# ./sem_utils/submit_slurm_jobs_sequential.sh stage01_structure_dcu/iter01/slurm model_update
