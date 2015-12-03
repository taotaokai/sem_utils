#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:?must provide control_file}
evid=${2:?must provide evid}

# source parameters in control_file
source ${control_file}

echo "Start plotting misfit: $(date)"

#====== setup event dir
event_dir=${iter_dir}/$evid

cd $event_dir

rm syn obs
ln -sf OUTPUT_forward syn
ln -sf $data_dir/$evid/disp obs

cd $event_dir

$wkdir/bin/plot_misfit.py obs syn misfit

echo "Done. $(date)"
