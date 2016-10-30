#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}
job_dep=${2:--1}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
source_dir=$wkdir/source

#------ source inversion
#work_flow=green
#work_flow=misfit
#work_flow=srcfrechet
#work_flow=dgreen
#work_flow=search

#------ structure inversion
#work_flow=syn
work_flow=misfit
#work_flow=kernel
#work_flow=hess
#work_flow=precond
#work_flow=perturb
#work_flow=search
#work_flow=hess_diag
#work_flow=hess_adj
#work_flow=hess_kernel

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do
  echo "====== $event_id"

  event_dir=$wkdir/$event_id

# #rm -rf $event_dir
# mkdir -p $event_dir/DATA

# # copy CMTSOLUTION file
# cmt_file=$(find -L $source_dir -name "${event_id}.cmt" | sort | tail -n1)
# echo $cmt_file
# if [ ! -f "$cmt_file" ]; then
#   echo "[ERROR] $cmt_file not found"
#   exit -1
# fi
# chmod u+w $event_dir/DATA/CMTSOLUTION.init
# cp $cmt_file $event_dir/DATA/CMTSOLUTION.init

# # copy STATIONS
# station_file=$data_dir/$event_id/data/STATIONS
# if [ ! -f "$station_file" ]; then
#   echo "[ERROR] $station_file not found"
#   exit -1
# fi
# #cp $station_file $event_dir/DATA/STATIONS
# #remove stations too close to the mesh western boundary
# awk '$1!~/#/{if($3<40 && $4>90) print $0; if($3>=40 && $4>=87) print $0;}' $station_file > $event_dir/DATA/STATIONS

# # copy misfit_par file
# cp $wkdir/misfit_par.py $event_dir/DATA/

# # create batch scripts
# #$utils_dir/make_source_iteration.sh $event_id
# $utils_dir/make_structure_iteration.sh $event_id

  $utils_dir/submit_slurm_jobs.sh $event_id ${work_flow} ${job_dep}

done