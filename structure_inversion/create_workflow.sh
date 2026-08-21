#!/bin/bash

#------ read command line args
control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}
workflow_file=${3:?[arg] workflow_file, e.g. workflow.json}
start_iter_num=${4:?[arg] start_iter_num, e.g. 0}
end_iter_num=${5:?[arg] end_iter_num, e.g. 5}

source ${control_file}
if [ $? -ne 0 ]
then
  echo "[ERROR] failed to source $control_file"
  exit -1
fi

cat <<EOF > ${workflow_file}
{
  "settings": {
    "max_retries": 5,
    "check_interval_seconds": 60,
    "max_concurrent_submissions": 100,
    "ban_duration_seconds": 3600,
    "max_node_strikes": 3,
    "log_file": "slurm.log",
    "checkpoint_path": "checkpoint.json"
  },
  "jobs": [
EOF

dep_job=

for iter in $(seq $start_iter_num $end_iter_num)
do
    
    # iteration directory
    SEM_iter_num=$(echo "$iter" | awk '{printf "%02d", $1}')
    SEM_iter_dir=${SEM_stage_dir}/iter${SEM_iter_num}
    SEM_iter_tag=s${SEM_stage_num}i${SEM_iter_num}

    mkdir -p $SEM_iter_dir/slurm
    cat<<EOF > $SEM_iter_dir/slurm/create_slurm_jobs.job
sed -i "0,/^SEM_iter_num/s/^SEM_iter_num=.*/SEM_iter_num=${SEM_iter_num}/" ${control_file}
$SEM_utils_dir/structure_inversion/proc_event_regEQ.sh ${control_file} ${event_list}
EOF

    job_name=${SEM_iter_tag}_create_jobs
    job_file=$SEM_iter_dir/slurm/create_slurm_jobs.job
    cat<<EOF >> $workflow_file
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": [$dep_job]
    }, 
EOF

    $SEM_utils_dir/structure_inversion/make_jobs_for_structure_iteration.sh \
        $SEM_iter_dir $SEM_iter_tag $event_list $workflow_file \"$job_name\"

    dep_job=\"${SEM_iter_tag}_model_reparam\"

done

cat <<EOF >> ${workflow_file}
  {}
  ]
}
EOF