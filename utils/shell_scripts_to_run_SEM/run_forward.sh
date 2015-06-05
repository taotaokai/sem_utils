#!/bin/bash

#prepare forward modelling for the current iteration
echo $(date)

# read Par_file
source Par_file

jobfile=forward.job

# check if mehser finished successfully
#mesher_status=$(awk -F"=" '/done/{print $2}' $iter_dir/mesher.status)
#if [ -z "$mesher_status" ] || [ "$mesher_status" -ne 1 ]
#then
#	echo "ERROR: Run run_mesher.sh first!"
#	exit
#fi

# prepare each event directory for forward simulation
#echo prepare event directory ...
#for evnm in $(cat $iter_dir/DATA/EVENTS)
#do
#	evnm_dir=$iter_dir/EVENTS/$evnm
#    echo $evnm_dir
#	cd $evnm_dir/DATABASES_MPI/
#	rm *
#	ln -s ../../../DATABASES_MPI/* ./ # link all mesh/topo files into local database directory
#done

# generate the job scrpit
echo make job file ...
cd $iter_dir
cat > $jobfile \
<<EOF
#$ -V                     		# Inherit the submission environment 
#$ -cwd                   		# Start job in submission directory
#$ -N forward.${iter}      		# Job Name
#$ -j y                   		# combine stderr & stdout into stdout  
#$ -o $iter_dir/forward.o\$JOB_ID       # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $num_proc    		# Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal              		# Queue name
#$ -l h_rt=$run_time_forward  # Run time (hh:mm:ss) - 1.5 hours
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
#
echo >> status
echo "Forward simulation" >> status
echo "begin: \$(date)" >> status
# change simulation type in DATA/
sed -i "/SIMULATION_TYPE/s/=.*/= 1/" DATA/Par_file
sed -i "/SAVE_FORWARD/s/=.*/= .true/" DATA/Par_file
#
rm OUTPUT_FILES
ln -s OUTPUT_forward OUTPUT_FILES
rm -rf OUTPUT_forward/*
cp $iter_dir/OUTPUT_mesher/addressing.txt OUTPUT_forward/
#
cd DATABASES_MPI/
find . -type l | xargs rm # remove all symlinks
ln -s $iter_dir/DATABASES_MPI/* ./ # link all mesh/topo files into local database directory
#
cd $evnm_dir/
ibrun $base_dir/SEM_bin/xspecfem3D
echo "end: \$(date)" >> status
#
echo >> status
echo "Make adjoint source" >> status
echo "begin: \$(date)" >> status
# set sac header
cd $evnm_dir/OUTPUT_forward
sac<<END
r *.sac
wh
q
END
# make adjoint sources
cd $evnm_dir
$USER_BIN/measure_adj_src.py DATA/STATIONS \
                             $event_dir/$evnm/DATA/PHASES_ADJOINT \
                             $event_dir/$evnm/DISPLACEMENT \
                             OUTPUT_forward \
                             ADJOINT_SOURCES \
                             $freqmin $freqmax > ADJOINT_SOURCES/adj.out
# sort adj.out with cc
sort -k14 -n ADJOINT_SOURCES/adj.out > tmp
mv tmp ADJOINT_SOURCES/adj.out
# generate STATION_ADJOINT
#choose only stations of cc >= cutoff(this value is now chosen arbitarily)
awk '\$14>c{printf "%5s %5s %8s %9s %8s   0\n",\$2,\$1,\$3,\$4,\$5}' c=$cc_cutoff \
    ADJOINT_SOURCES/adj.out >  DATA/STATIONS_ADJOINT
cp DATA/STATIONS_ADJOINT ADJOINT_SOURCES/

echo "end: \$(date)" >> status

EOF

done

# submit 
echo Now you can submmit job file by running
echo "> qsub $iter_dir/$jobfile"

#END
