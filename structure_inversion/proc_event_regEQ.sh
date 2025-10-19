#!/bin/bash
#set -x

# structure inversion for regional earthquake

# Work flows contain:
#--  kernel
#work_flow=forward,misfit,kernel
#-- line search
#work_flow=perturb,search
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

# create directories
mkdir -p $SEM_iter_dir

# check starting model exists
if [ ! -d "${SEM_starting_model_dir}" ]
then
  echo "[ERROR] ${SEM_starting_model_dir} does NOT exist!"
  exit -1
fi

# link initial model dir
initial_model_dir=${SEM_prev_iter_dir}/model_updated
if [ "$SEM_iter_num" -eq 0 ]
then
  initial_model_dir=${SEM_starting_model_dir}
fi
if [ -e $SEM_iter_dir/model_initial ] 
then
  rm -rf $SEM_iter_dir/model_initial
fi
ln -s $initial_model_dir $SEM_iter_dir/model_initial
# ln -sf -t $SEM_iter_dir/model_initial ${initial_model_dir}/*.bin 

# check Par_file exists
if [ ! -f "${SEM_config_dir}/DATA/Par_file" ]
then
  echo "[ERROR] ${SEM_config_dir}/DATA/Par_file does NOT exist!"
  exit -1
fi

# create slurm jobs for pre- and post-processing
$SEM_utils_dir/structure_inversion/make_slurm_jobs_for_pre_and_post_proc_regEQ.sh $control_file $event_list

# create slurm jobs for each events
# for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  echo "====== $event_id"

  # create event dir
  event_dir=$SEM_iter_dir/events/$event_id
  if [ ! -d "$event_dir" ]
  then
    mkdir -p $event_dir
  fi

  # copy Par_file
  mkdir -p $event_dir/DATA
  cp ${SEM_config_dir}/DATA/Par_file $event_dir/DATA/

  # copy latest updated CMTSOLUTION file
  # cmt_file=$(find -L $source_dir -path "*/iter??/events/${event_id}/misfit/CMTSOLUTION.updated" | sort | tail -n1)
  cmt_file=$(ls ${SEM_source_dir}/iter??/events/${event_id}/misfit/CMTSOLUTION.updated | sort | tail -n1)
  echo ------ source: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi
  [ -e $event_dir/DATA/CMTSOLUTION ] && chmod u+w $event_dir/DATA/CMTSOLUTION
  cp $cmt_file $event_dir/DATA/CMTSOLUTION

  # copy STATIONS
  station_file=$SEM_data_dir/$event_id/STATIONS
  if [ ! -f "$station_file" ]; then
    echo "[ERROR] $station_file not found"
    exit -1
  fi
  cp $station_file $event_dir/DATA/STATIONS

  # copy misfit_par file
  # cp $misfit_par_dir/${event_id}_misfit.yaml $event_dir/DATA/misfit.yaml
  cp $SEM_misfit_par_dir/misfit.yaml $event_dir/DATA/misfit.yaml

  # create slurm jobs
  $SEM_utils_dir/structure_inversion/make_slurm_jobs_for_event_regEQ.sh $control_file $event_id

done
