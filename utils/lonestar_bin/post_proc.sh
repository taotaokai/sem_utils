#!/bin/bash

# post-processing

#====== command line args
control_file=${1:?must provide control_file}
job=${2:?must prove job type(tar_forward,print_misfit)}
event_list=${3:?must provide evid list}

# log output
echo "# Start on $(date +%Y-%m-%dT%H:%M:%S)"
echo "# PWD: $HOSTNAME:$(pwd)"
echo "# COMMAND: $0 $@"
echo "# Parameter: control_file = $control_file"
echo "# Parameter: job = $job"
echo "# Parameter: event_list = $event_list"
echo

# load control parameters
source $control_file

#====== processes

#------ tar OUTPUT_forward
if [ "$job" == 'tar_forward' ]
then
    for evid in $(grep -v ^# $event_list)
    do
        echo ====== processing $evid
        event_dir=$iter_dir/$evid
        cd $event_dir
        ls OUTPUT_forward/* > /dev/null # touch all the files to avoid missing files in the tar
        tar -cf OUTPUT_forward.tar OUTPUT_forward
    done
fi

#------ print info
if [ "$job" == 'print_misfit' ]
then
    for evid in $(grep -v ^# $event_list)
    do
        echo "# $evid"
        event_dir=$iter_dir/$evid
        misfit_dir=$event_dir/misfit
        misfit_file=$misfit_dir/misfit.json
        tmp_file=$(mktemp)
        if [ -f ${misfit_file} ]
        then
            $wkdir/bin/print_misfit.py $misfit_file > $tmp_file
        else
            rm $tmp_file
            continue
        fi
        awk 'NF&&$1!~/^#/{dt=$7; cc0=$4; w=$11; 
            if(dt<0)dt=-dt; dts+=dt*w; ws+=w; cc0s+=cc0*w}
            END{printf "#%s Mean(dT/T) %g Mean(CC0) %g Sum(weight) %g\n",
            evid, dts/ws, cc0s/ws, ws}' evid=$evid  \
            $tmp_file > $misfit_dir/misfit.txt
        cat $tmp_file >> $misfit_dir/misfit.txt
        rm $tmp_file
    done
fi

# END
echo "#Finished on $(date +%Y-%m-%dT%H:%M:%S)"
