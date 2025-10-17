#!/bin/bash
#set -x

# structure inversion  for regional earthquake

# Work flows contain:
#--  kernel
#work_flow=syn,misfit,kernel
#work_flow=syn
#work_flow=misfit
#work_flow=kernel
#-- line search
#work_flow=perturb
#work_flow=search
#-- hessian-random model product
#work_flow=perturb_random,misfit_random,kernel_random
#work_flow=perturb_random
#work_flow=misfit_random
#work_flow=kernel_random

#------ read command line args
control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}

# load parameters in control_file
source $control_file

if [ ! -d "$stage_dir" ]
then
  mkdir -p $stage_dir
fi
if [ ! -d "$updated_model_dir" ]
then
  mkdir -p $updated_model_dir
fi

# initial model
echo ------ model: $(readlink -f ${initial_model_dir})
if [ ! -d "${initial_model_dir}" ]
then
  echo "[ERROR] ${initial_model_dir} does NOT exist!"
  exit -1
fi
if [ -d "$iter_dir/model_initial" ]
then
  rm -rf $iter_dir/model_initial
fi
ln -sf ${initial_model_dir} $iter_dir/model_initial

# sem Par_file
sem_par_file=${sem_config_dir}/DATA/Par_file
if [ ! -f "$sem_par_file" ]
then
  echo "[ERROR] $sem_par_file does NOT exist!"
  exit -1
fi

# mkdir -p $mesh_dir/DATA
# cp $sem_par_file $mesh_dir/DATA/Par_file
$sem_utils_dir/structure_inversion/make_slurm_jobs_for_pre_and_post_proc_regEQ.sh $control_file $event_list

#------ process each event
# for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  echo "====== $event_id"

  # create event dir
  event_dir=$iter_dir/events/$event_id
  sem_data_dir=$event_dir/DATA
  if [ -d "$sem_data_dir" ]
  then
    chmod u+w -R $event_dir/DATA
  fi
  mkdir -p $event_dir/DATA

  # copy Par_file
  cp $sem_par_file $event_dir/DATA/Par_file

  # copy CMTSOLUTION file
  # cmt_file=$(find -L $source_dir -path "*/iter??/events/${event_id}/misfit/CMTSOLUTION.updated" | sort | tail -n1)
  cmt_file=$(ls ${source_dir}/iter??/events/${event_id}/misfit/CMTSOLUTION.updated | sort | tail -n1)
  echo ------ source: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi
  [ -e $event_dir/DATA/CMTSOLUTION ] && chmod u+w $event_dir/DATA/CMTSOLUTION
  cp $cmt_file $event_dir/DATA/CMTSOLUTION

  # copy STATIONS
  station_file=$data_dir/$event_id/STATIONS
  if [ ! -f "$station_file" ]; then
    echo "[ERROR] $station_file not found"
    exit -1
  fi
  cp $station_file $event_dir/DATA/STATIONS

  # copy misfit_par file
  # cp $misfit_par_dir/${event_id}_misfit.yaml $event_dir/DATA/misfit.yaml
  cp $misfit_par_dir/misfit.yaml $event_dir/DATA/misfit.yaml

  # create batch scripts
  $sem_utils_dir/structure_inversion/make_slurm_jobs_for_event_regEQ.sh $control_file $event_id

done
