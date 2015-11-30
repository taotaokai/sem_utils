#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:?must provide control_file}
evid=${2:?must provide evid}

# source parameters in control_file
source ${control_file}

echo "Start measuring adjoint source: $(date)"

#====== setup event dir
event_dir=${iter_dir}/$evid

cd $event_dir

rm syn obs
ln -sf OUTPUT_forward syn
ln -sf $data_dir/$evid/disp obs

cd DATA
rm station.txt
ln -sf $data_dir/$evid/data/station.txt ./

cd ../
mkdir adj misfit

cd $event_dir
$wkdir/bin/measure_misfit.py DATA obs syn adj misfit $freqmin $freqmax

echo "Done. $(date)"
