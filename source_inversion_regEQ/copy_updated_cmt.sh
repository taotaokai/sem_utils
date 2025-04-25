#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

mkdir CMTSOLUTION_updated

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  # copy CMTSOLUTION
  cp $wkdir/$event_id/misfit/CMTSOLUTION.updated  $wkdir/CMTSOLUTION_updated/${event_id}.cmt

done