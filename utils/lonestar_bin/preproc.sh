#!/bin/bash

# Download files from jsg 

wkdir=$(pwd)

#====== command line args
iter=${1:?must provide iter#}
host=${2:-jsg20:~/NEChina}
event_list=${3:-evid}

# log output
echo "# Start on $(date +%Y-%m-%dT%H:%M:%S)"
echo "# PWD: $HOSTNAME:$(pwd)"
echo "# COMMAND: $0 $@"

#====== process
iter_dir=$(printf "iteration.%02d" $iter)

for evid in $(grep -v ^# $event_list); do

    echo
    echo ====== processing $evid
    echo

    event_dir=$wkdir/$iter_dir/$evid

    #----- remote copy misfit 
#   remote_file=$host/$iter_dir/$evid/misfit
#   rsync -auv $remote_file $event_dir/

    #----- remote copy adj 
    mkdir -p $event_dir/adj
    remote_file=$host/$iter_dir/$evid/adj
    rsync -auv $remote_file/*.adj $event_dir/adj
    
    #----- remote copy relocation
#   cd $event_dir
#   mkdir relocation; cd relocation
#   ln -sf $event_dir/misfit/relocate_1d_fix_depth.CMTSOLUTION CMTSOLUTION.new

done


# END
echo "#Finished on $(date +%Y-%m-%dT%H:%M:%S)"
