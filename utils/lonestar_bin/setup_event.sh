#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:?must provide control_file}
evid=${2:?must provide evid}

# source parameters in control_file
source ${control_file}

#====== create event dir
event_dir=${iter_dir}/$evid

if [ ! -d ${event_dir} ];then
    mkdir -p $event_dir
else
    echo "[ERROR] event_dir=$event_dir already exits!"
    exit -1
fi

#====== setup event dir
cd $event_dir
mkdir DATA DATABASES_MPI OUTPUT_forward

#---- setup DATA
cd $event_dir/DATA

#-- Par_file
cp -L $config_dir/DATA/Par_file ./

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
station_file=$data_dir/$evid/data/STATIONS
if [ -f $station_file ]
then
    cp -L $station_file ./
else
    echo "[ERROR] $station_file does NOT exist!"
    exit -1
fi

#-- CMTSOLUTION
raw_cmt_file=$data_dir/$evid/data/CMTSOLUTION
prev_cmt_file=$prev_iter_dir/$evid/misfit/CMTSOLUTION.reloc
if [ "$iter" -eq 0 ] || [ ! -f $prev_cmt_file ]
then
    # copy raw CMTSOLUTION
    if [ -f $raw_cmt_file ]
    then
        cp -L $raw_cmt_file CMTSOLUTION
    else
        echo "[ERROR] $raw_cmt_file does NOT exist!"
        exit -1
    fi
else
    cp -L $prev_cmt_file CMTSOLUTION
fi
# modify the header line to include time shift
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

#---- setup DATABASES_MPI
cd $event_dir/DATABASES_MPI

proc_file=$mesh_dir/DATABASES_MPI/proc000000_reg1_solver_data.bin
if [ ! -f "$proc_file" ]
then
    echo "[ERROR] Run mesher first. $proc_file does NOT exit!"
    exit -1
fi
ln -sf $mesh_dir/DATABASES_MPI/*.bin ./

#---- setup OUTPUT_FILES
cd $event_dir
ln -sf OUTPUT_forward OUTPUT_FILES
# copy addressing.txt from mesh_dir to OUTPUT_FILES
cp $mesh_dir/OUTPUT_FILES/addressing.txt $event_dir/OUTPUT_FILES

#---- backup parameter files
cd $event_dir/
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES

##====== make job script
#cat <<EOF > $event_dir/forward.job
##$ -V                              # Inherit the submission environment 
##$ -cwd                            # Start job in submission directory
##$ -N $evid.f$iter                 # Job Name
##$ -j y                            # combine stderr & stdout into stdout  
##$ -o $event_dir/forward.o         # Name of the output file (eg. myMPI.oJobID)
##$ -pe 12way $nproc_request        # Requests 12 cores/node, 24 cores total: 12way 24
##$ -q normal                       # Queue name
##$ -l h_rt=$run_time_forward       # Run time (hh:mm:ss) - 1.5 hours
##$ -M kai.tao@utexas.edu           # email 
##$ -m bea                          # email info: begin/end/abort
##$ -hold_jid -1                    # dependent job id
#
#export MY_NSLOTS=$nproc
#
#cd $event_dir/
#ibrun $build_dir/bin/xspecfem3D
#
#EOF
#
#echo "The forward modelling can be run as"
#echo ">> qsub $event_dir/forward.job"
