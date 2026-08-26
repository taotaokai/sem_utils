#!/bin/bash
#SBATCH -J make_workflow
#SBATCH -o make_workflow.job.o%j
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -p xahcnormal

# set -x

control_file=${1:?[arg]need control_file}
event_list=${2:?[arg] event_list, e.g. events.lst]}
workflow_file=${3:?[arg] workflow_file, e.g. workflow.json}
dep_jobs=${4:-}

source ${control_file}
if [ $? -ne 0 ]
then
  echo "[ERROR] failed to source $control_file"
  exit -1
fi

$SEM_utils_dir/structure_inversion/proc_event_regEQ.sh ${control_file} ${event_list}

SEM_iter_tag=s${SEM_stage_num}i${SEM_iter_num}

job_file=${SEM_iter_dir}/slurm/mesh.job
job_name=${SEM_iter_tag}_mesh
cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": [$dep_jobs]
    },
EOF

job_list=$(mktemp -p .) #submitted_job.lst
dep_job=${job_name}
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  event_dir=$SEM_iter_dir/events/$event_id
  prev_job=${dep_job}
  for job in forward misfit kernel reparam mask
  do
    job_file=${event_dir}/slurm/${job}.job
    job_name=${SEM_iter_tag}_${event_id}_${job}
    cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": ["$prev_job"]
    },
EOF
    prev_job=${job_name}
  done
  echo $job_name >> $job_list

  job_name=${SEM_iter_tag}_${event_id}_plot
  job_file=${event_dir}/slurm/plot.job
  cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": ["${SEM_iter_tag}_${event_id}_misfit"]
    },
EOF
done

dep_jobs=$(awk 'NR>1 {printf "\"%s\",",last } {last=$1} END{printf "\"%s\"",last}' $job_list)
for job in kernel_sum kernel_balance kernel_smooth dmodel model_perturb mesh_perturb
do
  job_file=${SEM_iter_dir}/slurm/${job}.job
  job_name=${SEM_iter_tag}_$job
  cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": [$dep_jobs]
    },
EOF
  dep_jobs=\"${job_name}\"
done

rm $job_list

job_list=$(mktemp -p .)
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  event_dir=$SEM_iter_dir/events/$event_id
  prev_job=$dep_jobs
  for job in perturb search
  do
    job_file=${event_dir}/slurm/${job}.job
    job_name=${SEM_iter_tag}_${event_id}_${job}
    cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": [$prev_job]
    },
EOF
    prev_job=\"${job_name}\"
  done
  echo $job_name >> $job_list
done

dep_jobs=$(awk 'NR>1 {printf "\"%s\",",last } {last=$1} END{printf "\"%s\"",last}' $job_list)
for job in model_update model_reparam
do
  job_file=${SEM_iter_dir}/slurm/${job}.job
  job_name=${SEM_iter_tag}_${job}
  cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": [$dep_jobs]
    },
EOF
  dep_jobs=\"${job_name}\"
done

rm $job_list

# cat <<EOF >> ${workflow_file}
#   {}
#   ]
# }
# EOF