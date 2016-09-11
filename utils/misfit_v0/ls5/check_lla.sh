#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  utils/ecef2lla.py CMTSOLUTION_initial/$event_id.cmt
  utils/ecef2lla.py CMTSOLUTION_updated/$event_id.cmt

done