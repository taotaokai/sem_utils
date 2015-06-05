#!/bin/bash

#prepare forward modelling for the current iteration
echo $(date)

# read Par_file
source Par_file

jobfile=adjoint.job

# check if mehser finished successfully
#mesher_status=$(awk -F"=" '/done/{print $2}' $iter_dir/mesher.status)
#if [ -z "$mesher_status" ] || [ "$mesher_status" -ne 1 ]
#then
#	echo "ERROR: Run run_mesher.sh first!"
#	exit
#fi

# check each event directory for adjoint simulation
echo check event directory ...
for evnm in $(cat $iter_dir/DATA/EVENTS)
do
	evnm_dir=$iter_dir/EVENTS/$evnm
    echo $evnm_dir
    # check if forward finished successfully
    #forward_status=$(awk -F"=" '/done/{print $2}' $evnm_dir/forward.status)
    #if [ -z "$forward_status" ] || [ "$forward_status" -ne 1 ]
    #then
   	#    echo ERROR: Run forward simutlation first!
    #    exit
    #fi

    # check if STATIONS_ADJOINT exists
    cd $evnm_dir/DATA
    if [ ! -f STATIONS_ADJOINT ]
    then
        echo "WARN: $evnm_dir/DATA/STATIONS_ADJOINT does not exist!"
    fi
done

# generate the job scrpit
echo make job file ...
cd $iter_dir
cat > $jobfile \
<<EOF
#$ -V                     		# Inherit the submission environment 
#$ -cwd                   		# Start job in submission directory
#$ -N adjoint.${iter}      		# Job Name
#$ -j y                   		# combine stderr & stdout into stdout  
#$ -o $iter_dir/adjoint.o\$JOB_ID       # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $num_proc    		# Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal              		# Queue name
#$ -l h_rt=$run_time_adjoint    # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu
#$ -m bea
#$ -hold_jid -1

# Run the MPI executable
EOF

for evnm in $(cat $iter_dir/DATA/EVENTS)
do
	evnm_dir=$iter_dir/EVENTS/$evnm
	cat >> $jobfile \
<<EOF
cd $evnm_dir/
echo >> status
echo "Adjoint simulation" >> status
echo "begin: \$(date)" >> status
# change simulation type in DATA/
sed -i "/SIMULATION_TYPE/s/=.*/= 3/" DATA/Par_file
sed -i "/SAVE_FORWARD/s/=.*/= .false/" DATA/Par_file
#
rm SEM
ln -s ADJOINT_SOURCES SEM
#
rm OUTPUT_FILES
ln -s OUTPUT_adjoint OUTPUT_FILES
rm -rf OUTPUT_adjoint/*
cp $iter_dir/OUTPUT_mesher/addressing.txt OUTPUT_adjoint/
#
cd $evnm_dir/
ibrun $base_dir/SEM_bin/xspecfem3D
echo "end: \$(date)" >> status

EOF

done

# submit 
echo Now you can submmit job file by running
echo "> qsub $iter_dir/$jobfile"

#END
