#!/bin/bash

# post-processing

#====== command line args
control_file=${1:?[arg] need control_file}
job=${2:?[arg] need job type(tar_forward,print_misfit)}
event_list=${3:?[arg] need event_list}

# log output
echo "# Start on $(date +%Y-%m-%dT%H:%M:%S)"
echo "# PWD: $HOSTNAME:$(pwd)"
echo "# COMMAND: $0 $@"
echo "# Parameter: control_file = $control_file"
echo "# Parameter: job = $job"
echo "# Parameter: event_list = $event_list"
echo

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)

if [ ! -f "$event_list" ]
then
    echo "[ERROR] invalid event_list: ", $event_list
    exit -1
fi
event_list=$(readlink -f $event_list)

# load control parameters
source $control_file

#====== processes

#------ tar OUTPUT_forward
if [ "$job" == 'tar_forward' ]
then
    for event_id in $(grep -v ^# $event_list)
    do
        echo "#====== processing $event_id"
        event_dir=$iter_dir/$event_id
        grep "End of the simulation" \
            $event_dir/OUTPUT_forward/output_solver.txt
        if [ "$?" -eq 0 ]
        then
            cd $event_dir
            ls OUTPUT_forward/* > /dev/null # touch all the files to avoid missing files in the tar
            tar -cf OUTPUT_forward.tar OUTPUT_forward
        else
            echo "[WARNING] ${event_id}: forward NOT existed/finished."
        fi
    done
fi

#------ print info
if [ "$job" == 'print_misfit' ]
then
    for event_id in $(grep -v ^# $event_list)
    do
        echo "#====== processing $event_id"
        event_dir=$iter_dir/$event_id
        misfit_dir=$event_dir/misfit
        misfit_file=$misfit_dir/misfit.json
        tmp_file=$(mktemp)
        if [ -f ${misfit_file} ]
        then
            $sem_utils/utils/print_misfit.py $misfit_file > $tmp_file
        else
            echo "[WARNING] ${event_id}: misfit NOT existed."
            rm $tmp_file
            continue
        fi
        awk 'NF&&$1!~/^#/{dt=$7; cc0=$4; w=$11; 
            if(dt<0)dt=-dt; dts+=dt*w; ws+=w; cc0s+=cc0*w}
            END{printf "#%s Mean(dT/T) %g Mean(CC0) %g Sum(weight) %g\n",
            event_id, dts/ws, cc0s/ws, ws}' event_id=$event_id  \
            $tmp_file > $misfit_dir/misfit.txt
        cat $tmp_file >> $misfit_dir/misfit.txt
        rm $tmp_file
    done
fi

# END
echo
echo "Finished on $(date +%Y-%m-%dT%H:%M:%S)"
echo
