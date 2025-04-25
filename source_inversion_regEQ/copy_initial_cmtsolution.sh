#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

source_dir=$wkdir/stage03.source

chmod u+w -R CMTSOLUTION_initial
rm -rf CMTSOLUTION_initial
mkdir CMTSOLUTION_initial

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  # copy CMTSOLUTION
  cmt_file=$(find -L $source_dir -name "${event_id}.cmt" | sort | tail -n1)
  if [ -z "$cmt_file" ]
  then
    echo "[WARN] CMTSOLUTION not found in $source_dir"
    continue
  fi
  cp $cmt_file CMTSOLUTION_initial/$event_id.cmt
  #sed -i "s/^event name:.*/event name:        $event_id/" CMTSOLUTION_initial/$event_id.cmt

done