#!/bin/bash

wkdir=$(pwd)

source_dir=${wkdir%iter??}

event_list=${1:?[arg]need event_list}

mkdir -p $wkdir/CMTSOLUTION_initial

echo $source_dir
for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do
  echo "====== $event_id"

  # copy CMTSOLUTION file
  cmt_file=$(find -L $source_dir -maxdepth 3 -path "*/iter??/CMTSOLUTION_updated/${event_id}.cmt" | sort | tail -n1)
  echo ------ use: $(readlink -f $cmt_file)
  if [ ! -f "$cmt_file" ]
  then
    echo "[ERROR] $cmt_file not found"
    exit -1
  fi
  cp $cmt_file $wkdir/CMTSOLUTION_initial/

done