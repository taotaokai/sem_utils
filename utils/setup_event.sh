#!/bin/bash

# setup event folder for SEM forward simulation

#====== command line args
control_file=${1:?[arg] need control_file}
event_id=${2:?[arg] need event_id}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
    exit 1
fi
control_file=$(readlink -f $control_file)

# load control parameters
source ${control_file}

#====== setup event dir

#------ create event dir
event_dir=${iter_dir}/$event_id

if [ ! -d "${event_dir}" ]
then
    mkdir -p $event_dir
else
    echo "[ERROR] event_dir=$event_dir already exits!"
    exit 1
fi

cd $event_dir
mkdir DATA DATABASES_MPI OUTPUT_forward

#------ setup DATA
cd $event_dir/DATA

#-- Par_file
cp -L $sem_config_dir/DATA/Par_file ./

sed -i "/SIMULATION_TYPE/s/=.*/= 1/" Par_file
if [ "$save_forward" == 'no' ]
then
    echo "* SAVE_FORWARD is set to .false."
    sed -i "/SAVE_FORWARD/s/=.*/= .false./" Par_file
    sed -i "/PARTIAL_PHYS_DISPERSION_ONLY/s/=.*/= .false./" Par_file
    sed -i "/UNDO_ATTENUATION/s/=.*/= .false./" Par_file
else
    echo "* SAVE_FORWARD is set to .true."
    sed -i "/SAVE_FORWARD/s/=.*/= .true./" Par_file
    sed -i "/PARTIAL_PHYS_DISPERSION_ONLY/s/=.*/= .false./" Par_file
    sed -i "/UNDO_ATTENUATION/s/=.*/= .true./" Par_file
fi

#-- STATIONS
station_file=$data_dir/$event_id/data/STATIONS
if [ -f $station_file ]
then
    cp -L $station_file ./
else
    echo "[ERROR] $station_file does NOT exist!"
    exit 1
fi

#-- CMTSOLUTION
raw_cmt_file=$data_dir/$event_id/data/CMTSOLUTION
prev_cmt_file=$prev_iter_dir/$event_id/misfit/CMTSOLUTION.reloc
if [ -f $prev_cmt_file ]
then
    cp -L $prev_cmt_file CMTSOLUTION
elif [ "${iter}" -le "${iter0}" ] && [ -f ${raw_cmt_file} ]
then
    cp -L $raw_cmt_file CMTSOLUTION
else
    echo "[ERROR] failed to set CMTSOLUTION for iteration: $iter"
    exit 1
fi
# modify header line to include time shift
tshift=$(awk '$0~/time shift/{print $3}' CMTSOLUTION)
otime=$(awk 'NR==1{printf "%s-%s-%s %s:%s:%s %s second", 
    $2,$3,$4,$5,$6,$7,a}' a=$tshift CMTSOLUTION)
info=$(awk 'NR==1{for(i=8;i<=NF;i++) printf "%s ",$i}' CMTSOLUTION)
ctime=$(date -u -d "$otime" +"%Y %m %d %H %M %S.%N")
mv CMTSOLUTION CMTSOLUTION.obs
echo "GCMT $ctime $info" > CMTSOLUTION.forward
awk 'NR>=2' CMTSOLUTION.obs >> CMTSOLUTION.forward
# set time shift to zero
echo "* time shift is set to zero."
sed -i "/time shift/s/:.*/:      0.0/" CMTSOLUTION.forward
# set half duration to zero
echo "* half duration is set to zero."
sed -i "/half duration/s/:.*/:   0.0/" CMTSOLUTION.forward
ln -sf CMTSOLUTION.forward CMTSOLUTION

#------ setup DATABASES_MPI
cd $event_dir/DATABASES_MPI

proc_file=$mesh_dir/DATABASES_MPI/proc000000_reg1_solver_data.bin
if [ ! -f "$proc_file" ]
then
    echo "[ERROR] Run mesher first. $proc_file does NOT exit!"
    exit 1
fi
ln -sf $mesh_dir/DATABASES_MPI/*.bin ./

#------ setup OUTPUT_FILES
cd $event_dir
ln -sf OUTPUT_forward OUTPUT_FILES
# copy addressing.txt from mesh_dir to OUTPUT_FILES
cp $mesh_dir/OUTPUT_FILES/addressing.txt $event_dir/OUTPUT_FILES

#------ backup parameter files
cd $event_dir/
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES

#END
