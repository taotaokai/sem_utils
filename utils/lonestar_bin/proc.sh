#!/bin/bash

# run different stages for iteration: new, mesh, event, qsub

wkdir=$(pwd)

#====== command line args
control_file=${1:?must provide control_file}
job=${2:?must prove job type(new,mesh,qsub_mesh,forward,qsub_forward,tar_forward,measure_adj)}
event_list=${3:?must provide evid list}

# log output
echo "# Start on $(date +%Y-%m-%dT%H:%M:%S)"
echo "# PWD: $HOSTNAME:$(pwd)"
echo "# COMMAND: $0 $@"
echo "# Parameter: control_file = $control_file"
echo "# Parameter: job = $job"
echo "# Parameter: event_list = $event_list"
echo

#====== process
# load control parameters
source $control_file

## setup new iteration
#if [ "$job" == 'new' ];then
#    bin/add_iteration.sh ${control_file}
#fi
#
## setup mesh 
#if [ "$job" == 'mesh' ];then
#    bin/setup_mesh.sh ${control_file}
#fi
#
## qsub mesh 
#if [ "$job" == 'qsub_mesh' ];then
#    if [ -f "$mesh_dir/mesher.o" ]; then
#        echo "[ERROR] SKIP: $mesh_dir/mesher.o exits."
#        continue
#    fi
#    qsub $mesh_dir/mesher.job   
#fi
#
## setup forward
#if [ "$job" == 'event' ];then
#    for evid in $(grep -v ^# $event_list); do
#        echo ====== processing $evid
#        bin/setup_event.sh $evid ${control_file}
#    done
#fi
#
## qsub forward 
#if [ "$job" == 'qsub_forward' ];then
#    for evid in $(grep -v ^# $event_list); do
#        echo ====== processing $evid
#
#        event_dir=$iter_dir/$evid
#        if [ -f "$event_dir/forward.o" ]; then
#            echo "[ERROR] SKIP: $event_dir/forward.o exits."
#            continue
#        fi
#        qsub $iter_dir/$evid/forward.job   
#    done
#fi
#
#
## setup adjoint
#if [ "$job" == 'adjoint' ];then
#    for evid in $(grep -v ^# $event_list); do
#        echo ====== processing $evid
#
#        event_dir=$iter_dir/$evid
#        if [ ! -f "$event_dir/forward.o" ]; then
#            echo "[ERROR] run forward first!"
#            continue
#        fi
#        if [ ! -d "$data_dir/$evid/disp" ]; then
#            echo "[ERROR] $data_dir/$evid/disp not exits!"
#            continue
#        fi
#        $wkdir/bin/measure_adjoint.sh $evid ${control_file}
#        $wkdir/bin/setup_adjoint.sh $evid ${control_file}
#    done
#fi
#
#
## qsub adjoint
#if [ "$job" == 'qsub_adjoint' ];then
#    for evid in $(grep -v ^# $event_list); do
#        echo ====== processing $evid
#        event_dir=$iter_dir/$evid
#        if [ -f "$event_dir/adjoint.o" ]; then
#            echo "[ERROR] $event_dir/adjoint.o exits, SKIP"
#            continue
#        fi
#        qsub $iter_dir/$evid/adjoint.job   
#    done
#fi
#
#
# tar OUTPUT_forward
if [ "$job" == 'tar_forward' ]
then
    for evid in $(grep -v ^# $event_list)
    do
        echo ====== processing $evid
        event_dir=$iter_dir/$evid
        #if [ ! -f "$event_dir/forward.o" ]; then
        #    echo "[ERROR] Run forward first. $event_dir/forward.o does NOT exit."
        #    continue
        #fi
        cd $event_dir
        ls OUTPUT_forward/* > /dev/null # touch all the files to avoid missing files in the tar
        tar -cf OUTPUT_forward.tar OUTPUT_forward
    done
fi
#
#
## kernel 
#if [ "$job" == 'kernel' ];then
#    for evid in $(grep -v ^# $event_list); do
#        echo ====== processing $evid
#        event_dir=$iter_dir/$evid
#        if [ ! -f "$event_dir/adjoint.o" ]; then
#            echo "[ERROR] run adjoint first!"
#            continue
#        fi
#        # make kernel mask
#        cd $event_dir
#        $sem_utils/bin/xsem_make_kernel_mask \
#            $nproc DATABASES_MPI OUTPUT_forward/source.vtk  DATABASES_MPI \
#            $source_mask_radius $stop_depth $pass_depth
#    done
#fi




# END
echo "#Finished on $(date +%Y-%m-%dT%H:%M:%S)"
