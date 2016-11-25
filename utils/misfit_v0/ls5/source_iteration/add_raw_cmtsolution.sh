#!/bin/bash

wkdir=$(pwd)

data_dir=$wkdir/events

utils=~/seiscode/sem_utils/utils/misfit_v0

event_list=${1:?[arg]need event_list}

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  $utils/change_cmt_to_ECEF.py $data_dir/$event_id/data/CMTSOLUTION CMTSOLUTION_initial/$event_id.cmt

done