#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

#!!! make sure these folders/links do exits and are what you want
mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
source_dir=$wkdir/CMTSOLUTION_initial

#====== make and submit jobs
#work_flow=green,misfit,srcfrechet,dgreen,search
#work_flow=green
#work_flow=misfit
#work_flow=srcfrechet
#work_flow=dgreen
work_flow=search

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  event_dir=$wkdir/$event_id

  if [[ "$work_flow" == "green"* ]]
  then
    chmod u+w -R $event_dir/DATA
    mkdir -p $event_dir/DATA

    # copy initial CMTSOLUTION file
    cmt_file=$source_dir/${event_id}.cmt
    echo ------ use: $(readlink -f $cmt_file)
    if [ ! -f "$cmt_file" ]; then
      echo "[ERROR] $cmt_file not found"
      exit -1
    fi
    cp $cmt_file $event_dir/DATA/CMTSOLUTION.init

    # copy misfit_par file
    cp $wkdir/misfit_par/${event_id}_misfit_par.py $event_dir/DATA/misfit_par.py

    # create batch scripts
    $utils_dir/make_source_iteration.sh $event_id
  fi

  $utils_dir/submit_slurm_jobs.sh $event_id ${work_flow}

done