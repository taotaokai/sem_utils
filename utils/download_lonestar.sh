#!/usr/bin/env bash

# Download simulation results and model from remote work computers (e.g. Lonestar)

#====== command line args
control_file=${1:?[args] need control_file}
job=${2:?[args] need jobs (comma separated string, e.g. model,kernel,event,all)}

#------ create associate array job_arr
declare -A job_arr
for tag in ${job//,/ }
do
    job_arr["$tag"]=1
done

#------ check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: " $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)

#------ need an third input argument when job includes "event|all"
if [[ ${job_arr["event"]} ]] || [[ ${job_arr["all"]} ]]
then
    event_list=${3:?[args] need event_id list}
    if [ ! -f "$event_list" ]
    then
        echo "[ERROR] invalid event_list: " $event_list
        exit -1
    fi
    event_list=$(readlink -f $event_list)
fi

#------ log output
echo "# Start on $(date +%Y-%m-%dT%H:%M:%S)"
echo "# PWD: $HOSTNAME:$(pwd)"
echo "# COMMAND: $0 $@"
echo "# Parameter: control_file = $control_file"
echo "# Parameter: job = ${job}"
echo "# Parameter: event_list = $event_list"
echo

# load control parameters
source $control_file
# remote path
remote_base_dir="lonestar:~/NEChina"
remote_iter_dir=$(printf "%s/iterations/iteration.%02d" $remote_base_dir $iter)

#====== process

#------ download model 
if [[ ${job_arr["model"]} ]] || [[ ${job_arr["all"]} ]]
then
    echo
    echo ====== download model
    echo
    mkdir -p ${model_dir}
    cd ${iter_dir}
    rsync -auv "$remote_iter_dir/model" ./
    #rsync -auv "$remote_iter_dir/mesh/DATA" .
    #rsync -auv "$remote_iter_dir/mesh/OUTPUT_FILES" .
fi

#------ download event
if [[ ${job_arr["event"]} ]] || [[ ${job_arr["all"]} ]]
then
    echo
    echo ====== download event
    echo
    
    for event_id in $(grep -v ^# $event_list)
    do
        echo 
        echo "EVENT: $event_id" 
        echo 
        event_dir=$iter_dir/$event_id
        if [ ! -d $event_dir ]
        then
            mkdir -p $event_dir
        fi
        cd $event_dir
        rsync -auvL "$remote_iter_dir/$event_id/OUTPUT_forward.tar" ./
        tar xf OUTPUT_forward.tar
        rm -rf misfit DATA
        rsync -auvL "$remote_iter_dir/$event_id/misfit" ./
        rsync -auvL "$remote_iter_dir/$event_id/DATA" ./
    done
fi

#------ download kernel
if [[ ${job_arr["kernel"]} ]] || [[ ${job_arr["all"]} ]]
then
    echo
    echo ====== download kernel 
    echo

    mkdir -p ${kernel_dir}
    cd ${iter_dir}
    rsync -auv $remote_iter_dir/kernel ./
fi

# END
echo "#Finished on $(date +%Y-%m-%dT%H:%M:%S)"
