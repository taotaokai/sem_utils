#!/bin/bash

event_list=${1:?[arg]event_list, e.g. "C202411201443A\nC202411201443A"}
stage_dir=${2:?[arg]stage_dir, e.g. stage00.source/ (find stage_dir/iter??/*/misfit/CMTSOLUTION.updated)}
out_dir=${3:?[arg]out_dir, e.g. CMTSOLUTION_updated/}

for event_id in $(awk 'NF&&$1!~/#/{print $1}' $event_list)
do
  echo "====== $event_id"

  # copy CMTSOLUTION file
  cmt_file=$(find -L $stage_dir -path "*/iter??/*/misfit/CMTSOLUTION.updated" | sort | tail -n1)
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi
  cp $cmt_file $out_dir/${event_id}.cmt

done
