#!/bin/bash

#prepare forward modelling for the current iteration
echo $(date)

# read Par_file
source Par_file

jobfile=update.job

# generate the job scrpit
echo make job file ...
cd $iter_dir
cat > $jobfile \
<<EOF
#$ -V                     		# Inherit the submission environment 
#$ -cwd                   		# Start job in submission directory
#$ -N update.${iter}      		# Job Name
#$ -j y                   		# combine stderr & stdout into stdout  
#$ -o $iter_dir/update.o\$JOB_ID       # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $num_proc    		# Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal              		# Queue name
#$ -l h_rt=$run_time_update     # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu
#$ -m bea
#$ -hold_jid -1

# Run the MPI executable
EOF

cat >> $jobfile \
<<EOF
cd $iter_dir/
echo >> status
echo "Update model" >> status
echo "begin: \$(date)" >> status

echo "- create source mask file" >> status
cd $iter_dir
for evnm in \$(cat $iter_dir/DATA/EVENTS)
do
  evnm_dir=$iter_dir/EVENTS/\$evnm
  cd \$evnm_dir
  rm DATABASES_MPI/*_mask_source.bin 
  $bin_dir/xmask_source $mask_width OUTPUT_forward/source.vtk $num_proc DATABASES_MPI DATABASES_MPI
done

cd $iter_dir
echo "done. \$(date)" >> status

echo "- prepare directory EVENT_KERNELS/" >> status
cd $iter_dir
rm -rf KERNEL_SUM EVENT_KERNELS
mkdir KERNEL_SUM EVENT_KERNELS
cd EVENT_KERNELS
for evnm in \$(cat $iter_dir/DATA/EVENTS)
do
  evnm_dir=$iter_dir/EVENTS/\$evnm
  ln -s \$evnm_dir/DATABASES_MPI \$evnm
done

cd $iter_dir
echo "done. \$(date)" >> status

echo "- sum event kernels" >> status
cd $iter_dir
cp DATA/EVENTS DATA/kernels_list.txt
#ibrun $base_dir/SEM_bin/xsum_preconditioned_kernels
ibrun $base_dir/SEM_bin/xsum_kernels

echo "done. \$(date)" >> status

echo "smooth event kernels" >> status
rm KERNEL_SUM/*_smooth.bin
#for kernm in alpha_kernel beta_kernel rho_kernel
#do
#  ibrun $base_dir/SEM_bin/xsmooth_sem $sigma_h $sigma_v \$kernm KERNEL_SUM DATABASES_MPI 
#done

#I don't want to smooth the kernel, so just symlink _kernel to _kernel_smooth
cd KERNEL_SUM/
ls *_kernel.bin | awk '{sub(/\.bin/,"",\$1);printf "ln -s %s.bin %s_smooth.bin\\n",\$1,\$1}' > ln-sh
bash ln-sh
cd $iter_dir/

echo "done. \$(date)" >> status

echo "- output new model files" >> status
ibrun $base_dir/SEM_bin/xadd_model_tiso_iso $step_fac KERNEL_SUM MODEL_new

echo "end: \$(date)" >> status
EOF

# submit 
echo Now you can submmit job file by running
echo "> qsub $iter_dir/$jobfile"

#END
