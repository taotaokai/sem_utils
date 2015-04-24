#!/bin/bash
#$ -V                     		# Inherit the submission environment 
#$ -cwd                   		# Start job in submission directory
#$ -N interp           		# Job Name
#$ -j y                   		# combine stderr & stdout into stdout  
#$ -o interp.o$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12           		# Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal              		# Queue name
#$ -l h_rt=03:00:00   # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu
#$ -m bea
#$ -hold_jid -1

# Run the MPI executable
mesh_dir=plot/mesh
out_dir=plot/model
model_str=vpv,vph,vsv,vsh,rho,eta

for profile in GUL TNC YCH
do
  for iter in 1 2 3 4 5 6
  do
    echo ----------------
    echo $(date)
    echo $profile $iter
    topo_dir=/scratch/03244/ktao/China_13s/SEM_iterations/$iter/DATABASES_MPI
    bin/interp_xyz $mesh_dir/$profile.xyz $model_str $topo_dir $topo_dir \
                   $out_dir/${profile}_iter${iter}.out
  done
done
