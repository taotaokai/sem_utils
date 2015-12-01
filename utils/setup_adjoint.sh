#!/bin/bash

# setup event folder for SEM adjoint simulation

#====== command line args
control_file=${1:?must provide control_file}
event_id=${2:?must provide event_id}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)

# load control parameters in control_file
source ${control_file}

echo
echo "setup_adjoint begins. [$(date)]"
echo

event_dir=${iter_dir}/$event_id

#====== measure adjoint source
echo
echo "#====== get adjoint source [$(date)]"
echo

cd $event_dir
rm syn obs
ln -sf OUTPUT_forward syn
ln -sf $data_dir/$event_id/disp obs
cd DATA
rm station.txt
ln -sf $data_dir/$event_id/data/station.txt ./

cd $event_dir
mkdir adj misfit
$base_dir/bin/measure_misfit.py DATA obs syn adj misfit $freqmin $freqmax

#====== setup files
echo
echo "#====== setup for adjoint simulation [$(date)]"
echo

#----- Par_file 
cd $event_dir/DATA
# adjoint simulation
sed -i "/SIMULATION_TYPE/s/=.*/= 3/" Par_file
sed -i "/SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/PARTIAL_PHYS_DISPERSION_ONLY/s/=.*/= .false./" Par_file
sed -i "/UNDO_ATTENUATION/s/=.*/= .true./" Par_file
# kernel Cijkl 
sed -i "/ANISOTROPIC_KL/s/=.*/= .true./" Par_file
sed -i "/SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
sed -i "/APPROXIMATE_HESS_KL/s/=.*/= .true./" Par_file
sed -i "/USE_FULL_TISO_MANTLE/s/=.*/= .true./" Par_file

#----- SEM
cd $event_dir
if [ ! -d adj ];then
    echo "[ERROR] $event_dir/adj does NOT exist!"
    exit -1
fi
ln -sf adj SEM

#----- STATIONS_ADJOINT
cd $event_dir/DATA
find ../adj -name "*XZ.adj" | sed "s/.*\///g; s/\..XZ.adj//; s/\./[ ]*/" > adj.grep
grep -f adj.grep STATIONS > STATIONS_ADJOINT

#----- OUTPUT_FILES
cd $event_dir
mkdir OUTPUT_adjoint
rm OUTPUT_FILES
ln -sf OUTPUT_adjoint OUTPUT_FILES
# copy addressing.txt from mesh_dir to OUTPUT_FILES
cp $mesh_dir/OUTPUT_FILES/addressing.txt $event_dir/OUTPUT_FILES

#------ copy input files
cd $event_dir/
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES


echo
echo "setup_adjoint ends. [$(date)]"
echo

#EOF
