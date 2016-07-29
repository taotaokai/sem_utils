#!/bin/bash

# copy date to jsg
wkdir=$(pwd)

jsg_host="jsg15:/home/u1/kt23545/NEChina/"

iter_dir=$(echo $wkdir | sed "s/.*\/\([^\/]*\/[^\/]*$\)/\1/")
#iter_dir="stage01.structure/iter01"

event_list=${1:?[arg]need event_list}

#====== <event_id>/
cat<<EOF > $wkdir/rsync_exclude.list
output_kernel/sac
output_hess/sac
DATABASES_MPI
misfit.pkl
adj_hess
adj_kernel
SEM
OUTPUT_FILES
source_mask
EOF

awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do
  echo "====== $event_id"
  event_dir=$wkdir/$event_id
  rsync -auvz $event_dir ${jsg_host}/${iter_dir} --exclude-from $wkdir/rsync_exclude.list
done

#====== kernel
for folder in kernel_sum
do
  echo "====== $folder"
  rsync -auvz $folder ${jsg_host}/${iter_dir}
done

#====== shell scripts
rsync -auv *.sh *.py *.txt *.pdf ${jsg_host}/${iter_dir}