#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:?[arg] need control_file}
event_id=${2:?[arg] need event_id}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: " $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)

# source parameters in control_file
source ${control_file}

echo "Start plotting misfit: $(date)"

#====== setup event dir
event_dir=${iter_dir}/$event_id

if [ ! -d "$event_dir" ]
then
    echo "[ERROR] invalid event_dir: " $event_dir
    exit -1
fi
cd $event_dir

rm syn obs
ln -sf OUTPUT_forward syn
ln -sf $data_dir/$event_id/disp obs

cd $event_dir

$sem_utils/utils/plot_misfit.py obs syn misfit

echo "Done. $(date)"
