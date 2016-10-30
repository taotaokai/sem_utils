#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}
job_dep=${2:--1}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
source_dir=$wkdir/source

#model_name=cosine_200km_010
model_name=cosine_200km_100

#------
#work_flow=perturb,hess_adj,hess_kernel
#work_flow=perturb
#work_flow=hess_adj
work_flow=hess_kernel

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do
  echo "====== $event_id"

  event_dir=$wkdir/$event_id

  #mv $event_dir/DATABASES_MPI $event_dir/forward_saved_frames
  #chmod a-w -R $event_dir/forward_saved_frames

  # clean previous jobs
  rm $event_dir/slurm/${work_flow}.job*

  $utils_dir/make_structure_hessian.sh $event_id ${model_name}

  $utils_dir/submit_slurm_jobs.sh $event_id ${work_flow} ${job_dep}

done