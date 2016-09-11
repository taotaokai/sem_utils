#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
source_dir=$wkdir/source

#====== make and submit jobs

#------ source inversion
#work_flow=green
#work_flow=misfit
#work_flow=srcfrechet
#work_flow=dgreen
#work_flow=search

#------ structure inversion
work_flow=syn
#work_flow=misfit
#work_flow=kernel
#work_flow=hess
#work_flow=precond
#work_flow=perturb
#work_flow=search

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do
  echo "====== $event_id"

  event_dir=$wkdir/$event_id

  # copy initial CMTSOLUTION file
  cmt_file=$source_dir/${event_id}.cmt
  if [ ! -f "$cmt_file" ]; then
    echo "[ERROR] CMTSOLUTION not found"
    exit -1
  fi
  mkdir -p $event_dir/DATA
  cp $cmt_file $event_dir/DATA/CMTSOLUTION.init

  # copy misfit_par file
  cp $wkdir/misfit_par.py $event_dir/DATA/

  # create batch scripts
  $utils_dir/make_source_iteration.sh $event_id

  $utils_dir/submit_slurm_jobs.sh $event_id ${work_flow}

done