#!/bin/bash

# check job running status based on slurm output logs

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}
work_type=${2:?[arg]need work_type e.g. syn kernel hess}

#====== check error in stage syn 
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  #echo "====== $event_id"
  event_dir=$wkdir/$event_id

  #------ check if log file exists
  ls $event_dir/slurm/${work_type}.job.o* >/dev/null 2>/dev/null
  stat=$?
  # if log file not found
  if [ $stat -ne 0 ]; then
    echo $event_id
    continue 
  fi

  #------ check if any error in the last log file
  slurm_log=$(ls $event_dir/slurm/${work_type}.job.o* 2>/dev/null | sort | tail -n1)
  grep -i -e "error" -e "cancelled" -e "bad" $slurm_log > /dev/null
  stat=$?
  # if error found
  if [ $stat -eq 0 ]; then
    echo $event_id
    continue 
  fi

done

cat<<EOF
#Use this command to comment out those bad events from to_run.txt
#awk '{printf "sed -i \"/%s/s/^/#>/\" to_run.txt\n", \$1}' tmp
EOF