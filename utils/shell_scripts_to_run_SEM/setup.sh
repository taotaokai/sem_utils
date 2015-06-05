#!/bin/bash

#20150131 created
# setup folders for one iteration 

# get parameters
source Par_file

# create first-level folders
mkdir $iter_dir
cd $iter_dir
mkdir DATA DATABASES_MPI EVENTS GRADIENT UPDATE_DIRECTION OUTPUT_mesher
if [ $iter -gt 0 ]
then
  mkdir MODEL DIFF_MODEL DIFF_GRADIENT
else
  ln -s $base_dir/INITIAL_MODEL MODEL
fi


# copy/link files to DATA/
cd $iter_dir/DATA/
cp $config_dir/DATA/Par_file Par_file
ln -s $iter_dir/MODEL GLL
ln -s $config_dir/DATA/crust2.0
ln -s $config_dir/DATA/s362ani
ln -s $config_dir/DATA/QRFSI12
ln -s $config_dir/DATA/topo_bathy


# setup EVENTS/
cd $iter_dir/EVENTS/
cp $base_dir/EVENTS/LIST $iter_dir/DATA/EVENTS
for evnm in $(cat $iter_dir/DATA/EVENTS)
do
    evnm_dir=$iter_dir/EVENTS/$evnm

    # create folders
    mkdir $evnm_dir
    cd $evnm_dir
    mkdir DATA DATABASES_MPI ADJOINT_SOURCES CMT_new OUTPUT_adjoint OUTPUT_forward 

    # copy Par_file
    cp $iter_dir/DATA/Par_file $evnm_dir/DATA/

    # copy CMTSOLUTION
    cd $event_dir/$evnm/DATA
    iter_CMT=$(ls CMTSOLUTION.? | awk -F. '$2<a&&b<$2{b=$2}END{print b}' a=$iter b=0)
    if [ -z "$iter_CMT" ]
    then
        echo "ERROR: no available CMTSOLUTION in $event_dir/$evnm !"
        exit
    fi
    cp CMTSOLUTION.${iter_CMT} $evnm_dir/DATA/CMTSOLUTION

    # copy STATIONS
    if [ -f $event_dir/$evnm/DATA/STATIONS ]
    then
        cp $event_dir/$evnm/DATA/STATIONS $evnm_dir/DATA/
    else
        echo "ERROR: $event_dir/$evnm/DATA/STATIONS does not exist!"
        exit
    fi
done

# mesher should not depend on CMTSOLUTION,
# however, meshfem3D/initialize_mesher.f90 calls shared/read_compute_parameters.f90,
# which calls shared/count_number_of_sources.f90 that tries to open DATA/CMTSOLUTION.
# So to circumvent this problem, just put any valid CMTSOULTION under iter_dir/DATA/
cp $event_dir/$evnm/DATA/CMTSOLUTION.${iter_CMT} $iter_dir/DATA/CMTSOLUTION

#END
