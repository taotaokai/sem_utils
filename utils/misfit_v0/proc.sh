#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
source_dir=$wkdir/source

#jsg_host="jsg19:~/NEChina/source_inversion/"

#====== make and submit jobs
work_flow=syn,misfit,kernel,hess
#work_flow=misfit,kernel,hess

awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo "====== $event_id"

  $utils_dir/make_structure_iteration.sh $event_id

  # copy CMTSOLUTION
  cmt_file=$(ls $source_dir/$event_id/DATA/CMTSOLUTION.iter?? | sort | tail -n1)
  if [ $? -ne 0 ]; then
    echo "[ERROR] CMTSOLUTION not found"
    exit -1
  fi
  mkdir -p $wkdir/$event_id/DATA/
  cp $cmt_file $wkdir/$event_id/DATA/CMTSOLUTION
  sed -i "s/^event name:.*/event name:        $event_id/" $wkdir/$event_id/DATA/CMTSOLUTION

  # copy misfit_par file
  cp $utils_dir/misfit_par.py $wkdir/$event_id/DATA/

  $utils_dir/submit_structure_iteration.sh $event_id ${work_flow}

done