#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}
job_dep=${2:--1}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
source_dir=$wkdir/source

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do
  echo "#====== $event_id"

  event_dir=$wkdir/$event_id

  nfile=$(/bin/ls $event_dir/forward_saved_frames/proc000*_save_frame_at0000* | wc -l)

  if [ $nfile -ne 10752 ]
  then
    echo "$event_id $nfile"
  fi

done