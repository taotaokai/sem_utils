#!/bin/bash

# run mesher for the current iteration
echo $(date)

# read Par_file
source Par_file

jobfile=mesher.job

# generate job script
cd $iter_dir
cat > $jobfile \
<<EOF
#$ -V                     		# Inherit the submission environment 
#$ -cwd                   		# Start job in submission directory
#$ -N mesher.${iter}      		# Job Name
#$ -j y                   		# combine stderr & stdout into stdout  
#$ -o $iter_dir/mesher.o\$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $num_proc    		# Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal              		# Queue name
#$ -l h_rt=$run_time_mesher   # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu
#$ -m bea
#$ -hold_jid -1

# Run the MPI executable
cd $iter_dir

echo >> status
echo "Mesher" >> status
echo "begin: \$(date)" >> status

# create links
rm OUTPUT_FILES
ln -s OUTPUT_mesher OUTPUT_FILES

ibrun $base_dir/SEM_bin/xmeshfem3D

echo "end: \$(date)" >> status

EOF

# submit 
echo Now you can submit the job by running
echo "> qsub $iter_dir/$jobfile"

#END
