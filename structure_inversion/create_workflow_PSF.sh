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
initial_model_dir=${4:?[arg] initial_model_dir, e.g. mesh_initial/GLL_initial}
perturbed_model_dir=${5?:[arg] perturbed_model_dir, e.g. mesh_initial/GLL_initial_psf_alpha}
data_type=${6?:[arg]data_type, e.g. DISP, VEL}
dep_jobs=${7:-}

source ${control_file}
if [ $? -ne 0 ]
then
  echo "[ERROR] failed to source $control_file"
  exit -1
fi

if [ ! -d "${initial_model_dir}" ]
then
  echo "[ERROR] initial_model_dir(${initial_model_dir}) does NOT exist"
  exit -1
fi
initial_model_dir=$(readlink -f $initial_model_dir)

if [ ! -d "${perturbed_model_dir}" ]
then
  echo "[ERROR] perturbed_model_dir(${perturbed_model_dir}) does NOT exist"
  exit -1
fi
perturbed_model_dir=$(readlink -f $perturbed_model_dir)

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

#1. get synthetic seismogram from perturbed model 

sed -i "0,/^SEM_iter_num/s/^SEM_iter_num=.*/SEM_iter_num=0/" ${control_file}
$SEM_utils_dir/structure_inversion/proc_event_regEQ.sh ${control_file} ${event_list}

source ${control_file}

rm $SEM_iter_dir/model_initial
ln -sf ${perturbed_model_dir} $SEM_iter_dir/model_initial

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

  # run forward simulation
  job_file=${event_dir}/slurm/forward.job
  job_name=${SEM_iter_tag}_${event_id}_forward
  cat <<EOF >> ${workflow_file}
  {
    "name": "${job_name}",
    "script": "$job_file",
    "dependencies": ["$prev_job"]
  },
EOF

  # convert sac fils to .h5 file
  job_name=${SEM_iter_tag}_${event_id}_sac2h5
  job_file=${event_dir}/slurm/sac2h5.job
  cat <<EOF > ${job_file}
#!/bin/bash
#SBATCH -J ${job_name}
#SBATCH -o ${job_file}.o%j
#SBATCH -n 1
#SBATCH -p ${SLURM_partition_cpu}
#SBATCH -t 00:01:00

sac_dir=${event_dir}/output_forward/sac
station_file=${event_dir}/DATA/STATIONS
cmt_file=${SEM_data_dir}/${event_id}/CMTSOLUTION
out_h5=${event_dir}/data.h5
channel_file=${event_dir}/channel.txt

python $SEM_utils_dir/data_process/sem_sac2h5.py \
  \$sac_dir \$station_file \$out_h5 \
  --cmt_file \$cmt_file \
  --pad_before 600 \
  --data_type ${data_type}

awk '
BEGIN {
  print "#Network|Station|Location|Channel|Latitude|Longitude|Elevation|Depth|Azimuth|Dip|SensorDescription|Scale|ScaleFreq|ScaleUnits|SampleRate|StartTime|EndTime"
}
{
  split(\$2,parts,".");
  sta=parts[1]; loc=parts[2];
  printf "%s|%s|%s|BXZ|%f|%f|%f|%f|0|-90||||||1999-01-01|2030-01-01\n",\$1,sta,loc,\$3,\$4,\$5,\$6;
  printf "%s|%s|%s|BXE|%f|%f|%f|%f|90|0||||||1999-01-01|2030-01-01\n",\$1,sta,loc,\$3,\$4,\$5,\$6;
  printf "%s|%s|%s|BXN|%f|%f|%f|%f|0|0||||||1999-01-01|2030-01-01\n",\$1,sta,loc,\$3,\$4,\$5,\$6;
}' \$station_file > \$channel_file

cp \$station_file ${event_dir}/STATIONS 
EOF

  cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": ["${SEM_iter_tag}_${event_id}_forward"]
    },
EOF
  echo $job_name >> $job_list

done

#2. calculate kernel from initial model with perturbed model synthetics as data

SEM_data_dir=${SEM_iter_dir}/events

sed -i "0,/^SEM_iter_num/s/^SEM_iter_num=.*/SEM_iter_num=1/" ${control_file}
sed -i "0,/^SEM_data_dir/s/^SEM_data_dir=.*/SEM_data_dir=${SEM_data_dir}/" ${control_file}
$SEM_utils_dir/structure_inversion/proc_event_regEQ.sh ${control_file} ${event_list}

source ${control_file}

rm $SEM_iter_dir/model_initial
ln -sf ${initial_model_dir} $SEM_iter_dir/model_initial

dep_jobs=$(awk 'NR>1 {printf "\"%s\",",last } {last=$1} END{printf "\"%s\"",last}' $job_list)

job_file=${SEM_iter_dir}/slurm/mesh.job
job_name=${SEM_iter_tag}_mesh
cat <<EOF >> ${workflow_file}
    {
      "name": "${job_name}",
      "script": "$job_file",
      "dependencies": [$dep_jobs]
    },
EOF

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
for job in kernel_sum # kernel_balance kernel_smooth dmodel model_perturb mesh_perturb
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

cat <<EOF >> ${workflow_file}
  {}
  ]
}
EOF

rm $job_list