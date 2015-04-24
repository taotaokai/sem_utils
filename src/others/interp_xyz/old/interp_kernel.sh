#!/bin/bash
#$ -V                     		# Inherit the submission environment 
#$ -cwd                   		# Start job in submission directory
#$ -N interp_ker           		# Job Name
#$ -j y                   		# combine stderr & stdout into stdout  
#$ -o interp_ker.o$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12           		# Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal              		# Queue name
#$ -l h_rt=03:00:00   # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu
#$ -m bea
#$ -hold_jid -1

# Run the MPI executable
mesh_dir=plot/mesh
model_str=alpha_kernel,beta_kernel,rho_kernel

out_dir=plot/event_kernel
for profile in GUL TNC YCH
do
  for iter in 1 2 3 4 5 6
  do
    echo ----------------
    echo $(date)
    echo $profile $iter
    #topo_dir=/scratch/03244/ktao/China_13s/SEM_iterations/$iter/DATABASES_MPI
    topo_dir=/scratch/03244/ktao/China_13s/SEM_iterations/$iter/EVENT_KERNELS/20030727_0625
    bin/interp_xyz $mesh_dir/$profile.xyz $model_str $topo_dir $topo_dir \
                   $out_dir/${profile}_iter${iter}.out
  done
done


out_dir=plot/kernel_sum
for profile in GUL TNC YCH
do
  for iter in 1 2 3 4 5 6
  do
    echo ----------------
    echo $(date)
    echo $profile $iter
    topo_dir=/scratch/03244/ktao/China_13s/SEM_iterations/$iter/DATABASES_MPI
    model_dir=/scratch/03244/ktao/China_13s/SEM_iterations/$iter/KERNEL_SUM
    bin/interp_xyz $mesh_dir/$profile.xyz $model_str $model_dir $topo_dir \
                   $out_dir/${profile}_iter${iter}.out
  done
done
