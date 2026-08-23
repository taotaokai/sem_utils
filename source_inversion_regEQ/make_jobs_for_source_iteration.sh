#!/bin/bash
#SBATCH -J make_workflow
#SBATCH -o make_workflow.job.o%j
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -p xahcnormal

# set -x

# SEM_iter_dir=${1:?[arg] SEM_iter_dir, e.g. stage00_source/iter00}
# SEM_iter_tag=${2:?[arg] SEM_iter_tag, e.g. s00i01, job name prefix}
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

$SEM_utils_dir/source_inversion_regEQ/proc_event_regEQ.sh ${control_file} ${event_list}

SEM_iter_tag=s${SEM_stage_num}i${SEM_iter_num}

job_list=$(mktemp -p .) #submitted_job.lst

for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  event_dir=$SEM_iter_dir/events/$event_id
  prev_job=${dep_jobs}
  for job in green misfit srcfrechet dgreen search
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
job_file=${SEM_iter_dir}/slurm/sync.job
job_name=${SEM_iter_tag}_sync

cat <<EOF >> ${job_file}
#!/bin/bash
#SBATCH -J ${job_name}
#SBATCH -o ${job_file}.o%j
#sbatch -n 1
#sbatch -p ${SLURM_partition_cpu}
#sbatch -t 00:01:00
echo ${job_name}
EOF

cat <<EOF >> ${workflow_file}
  {
    "name": "${job_name}",
    "script": "$job_file",
    "dependencies": [$dep_jobs]
  },
EOF

rm $job_list

# cat <<EOF >> ${workflow_file}
#   {}
#   ]
# }
# EOF
