#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:?must provide control_file}
evid=${2:?must provide evid}

# source parameters in control_file
source ${control_file}

#====== setup event dir
event_dir=${iter_dir}/$evid

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

#====== make job script
#cat <<EOF > $event_dir/adjoint.job
##$ -V                              # Inherit the submission environment 
##$ -cwd                            # Start job in submission directory
##$ -N $evid.a$iter                 # Job Name
##$ -j y                            # combine stderr & stdout into stdout  
##$ -o $event_dir/adjoint.o$JOB_ID  # Name of the output file (eg. myMPI.oJobID)
##$ -pe 12way $nproc_request        # Requests 12 cores/node, 24 cores total: 12way 24
##$ -q normal                       # Queue name
##$ -l h_rt=$run_time_adjoint       # Run time (hh:mm:ss) - 1.5 hours
##$ -M kai.tao@utexas.edu           # email 
##$ -m bea                          # email info: begin/end/abort
##$ -hold_jid -1                    # dependent job id
#
#export MY_NSLOTS=$nproc

cd $event_dir/
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES

#ibrun $build_dir/bin/xspecfem3D

#EOF
#
##
#echo "The adjoint simulation can be run as"
#echo "qsub $event_dir/adjoint.job"
