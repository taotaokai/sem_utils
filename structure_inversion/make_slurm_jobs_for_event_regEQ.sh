#!/bin/bash

# Make jobs files for slurm
# structure inversion

control_file=${1:?[arg]need control_file}
event_id=${2:?[arg]need event_id}

# source control_file
if [ ! -f "$control_file" ]
then
  echo "[ERROR] $control_file NOT found!"
  exit -1
fi
source $control_file

#====== define variables

# directories
event_dir=$SEM_iter_dir/events/$event_id
mesh_dir=$SEM_iter_dir/mesh
misfit_dir=$event_dir/misfit
figure_dir=$misfit_dir/figure
slurm_dir=$event_dir/slurm
# job scripts for slurm
mkdir -p $slurm_dir
forward_job=$slurm_dir/forward.job
misfit_job=$slurm_dir/misfit.job
kernel_job=$slurm_dir/kernel.job
# process_job=$slurm_dir/process.job
reparam_job=$slurm_dir/reparam.job
mask_job=$slurm_dir/mask.job
threshold_job=$slurm_dir/threshold.job
plot_job=$slurm_dir/plot.job
#hess_job=$slurm_dir/hess.job
#precond_job=$slurm_dir/precond.job
perturb_job=$slurm_dir/perturb.job
search_job=$slurm_dir/search.job
##hess_diag_job=$slurm_dir/hess_diag.job
##hess_model_product_job=$slurm_dir/hess_model_product.job
##hess_kernel_job=$slurm_dir/hess_kernel.job
## hessian-random model product
#hess_syn_job=$slurm_dir/hess_syn.job
#hess_misfit_job=$slurm_dir/hess_misfit.job
#hess_kernel_job=$slurm_dir/hess_kernel.job

# database file
#mkdir -p $misfit_dir
db_file=$misfit_dir/misfit.h5

# station file
station_file=$event_dir/DATA/STATIONS
if [ ! -f "$station_file" ]
then
  echo "[ERROR] $station_file does NOT exist!"
  exit -1
fi

# cmt file
cmt_file=$event_dir/DATA/CMTSOLUTION
if [ ! -f "$cmt_file" ]
then
  echo "[ERROR] $cmt_file does NOT exist!"
  exit -1
fi

# misfit par file
misfit_par_file=$event_dir/DATA/misfit.yaml
if [ ! -f "$misfit_par_file" ]
then
  echo "[ERROR] $misfit_par_file does NOT exist!"
  exit -1
fi

# SEM Par_file
sem_par_file=${event_dir}/DATA/Par_file
if [ ! -f "$sem_par_file" ]
then
  echo "[ERROR] $sem_par_file does NOT exist!"
  exit -1
fi
# set USER_T0 in Par_file to at least 3 * tau
tau=$(grep "tau" ${cmt_file} | awk -F: '{printf "%f", $2}')
min_user_t0=$(echo "3 * $tau + 1" | bc -l | awk '{printf "%d", $1}')
sed -i "/^T0/s/=.*/= $min_user_t0/" ${sem_par_file}

#====== forward simulation
cat <<EOF > $forward_job
#!/bin/bash
#SBATCH -J ${event_id}.forward
#SBATCH -o $forward_job.o%j
#SBATCH ${SLURM_args_forward}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

out_dir=output_forward

mkdir -p $event_dir/DATA
cd $event_dir/DATA

sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .true./" Par_file

# sed -i "/^USE_ECEF_COORDINATE/s/=.*/= .true./" Par_file
# for regional earthquake data
# sed -i "/^USE_FORCE_POINT_SOURCE/s/=.*/= .false./" Par_file

rm -rf $event_dir/DATABASES_MPI
mkdir $event_dir/DATABASES_MPI
ln -s $mesh_dir/DATABASES_MPI/*.bin $event_dir/DATABASES_MPI

cd $event_dir
rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

${SLURM_mpiexec} $SEM_build_dir/bin/xspecfem3D

# check if simulation finished successfully
if grep -qF "End of the simulation" \$out_dir/output_solver.txt; then
  echo "====== forward job SUCCESS: $event_id"
else
  echo "====== forward job FAILED: $event_id"
  exit 1
fi

mkdir $event_dir/\$out_dir/sac
mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

chmod a+w -R $event_dir/forward_saved_frames
rm -rf $event_dir/forward_saved_frames
mv $event_dir/DATABASES_MPI $event_dir/forward_saved_frames
chmod a-w -R $event_dir/forward_saved_frames

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

#====== misfit
cat <<EOF > $misfit_job
#!/bin/bash
#SBATCH -J ${event_id}.misfit
#SBATCH -o ${misfit_job}.o%j
#SBATCH ${SLURM_args_misfit}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

export OPENBLAS_NUM_THREADS=1

cd $event_dir

chmod u+w -R $misfit_dir
rm -rf $misfit_dir
mkdir -p $misfit_dir

chmod u+w -R $event_dir/SEM
rm -rf $event_dir/SEM
mkdir -p $event_dir/SEM

$SEM_python_exec $SEM_utils_dir/misfit/measure_adj.py \\
  $db_file \\
  $misfit_par_file \\
  $sem_par_file \\
  $cmt_file \\
  $SEM_data_dir/$event_id/channel.txt \\
  $SEM_data_dir/$event_id/data.h5 \\
  $event_dir/output_forward/sac \\
  $event_dir/SEM \\
  --nproc=\${SLURM_NTASKS} \\
  --cmt_in_ECEF \\
  --window_yaml $event_dir/DATA/window.yaml

if [ \$? -ne 0 ]
then
  echo "measure_adj.py failed!"
  exit 1
fi

# make STATIONS_ADJOINT
cd $event_dir/SEM
ls *Z.adj | sed 's/..Z\.adj$//' |\
  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
grep -f $event_dir/SEM/grep_pattern $event_dir/DATA/STATIONS \
  > $event_dir/SEM/STATIONS_ADJOINT

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

#====== kernel simulation
cat <<EOF > $kernel_job
#!/bin/bash
#SBATCH -J ${event_id}.kernel
#SBATCH -o $kernel_job.o%j
#SBATCH ${SLURM_args_kernel}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

out_dir=output_kernel

cd $event_dir/DATA
chmod u+w Par_file
sed -i "/^SIMULATION_TYPE/s/=.*/= 3/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^ANISOTROPIC_KL/s/=.*/= .true./" Par_file
sed -i "/^SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
sed -i "/^APPROXIMATE_HESS_KL/s/=.*/= .false./" Par_file

cp -f $event_dir/SEM/STATIONS_ADJOINT $event_dir/DATA/
# rm -rf $event_dir/SEM
# ln -s $event_dir/adj_kernel $event_dir/SEM

cd $event_dir

rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

rm -rf $event_dir/DATABASES_MPI
mkdir $event_dir/DATABASES_MPI
ln -s $event_dir/forward_saved_frames/*.bin $event_dir/DATABASES_MPI

# cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

cd $event_dir
${SLURM_mpiexec} $SEM_build_dir/bin/xspecfem3D

# check if simulation finished successfully
if grep -qF "End of the simulation" \$out_dir/output_solver.txt; then
  echo "====== kernel job SUCCESS: $event_id"
else
  echo "====== kernel job FAILED: $event_id"
  exit 1
fi

mkdir $event_dir/\$out_dir/sac
mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

mkdir $event_dir/\$out_dir/kernel
mv $event_dir/DATABASES_MPI/*_kernel.bin $event_dir/\$out_dir/kernel/
# mv $event_dir/DATABASES_MPI/*reg1_cijkl_kernel.bin $event_dir/\$out_dir/kernel/
# mv $event_dir/DATABASES_MPI/*reg1_rho_kernel.bin $event_dir/\$out_dir/kernel/

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo
EOF


#====== reduce cijkl,rho kernel to VTI parameters
cat <<EOF > $reparam_job
#!/bin/bash
#SBATCH -J ${event_id}.reparam
#SBATCH -o $reparam_job.o%j
#SBATCH ${SLURM_args_kernel_process}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

echo ====== check if simulation finished successfully

kernel_dir=${event_dir}/output_kernel/kernel

n=\$(ls -1 \$kernel_dir/proc*_cijkl_kernel.bin | wc -l)
if [ "\$n" -ne "${SEM_nproc_total}" ]; then
  echo "number of \$kernel_dir/proc*_cijkl_kernel.bin files not equal to ${SEM_nproc_total}"
  echo "====== kernel process job FAILED: $event_id"
  exit 1
fi

n=\$(ls -1 \$kernel_dir/proc*_rho_kernel.bin | wc -l)
if [ "\$n" -ne "${SEM_nproc_total}" ]; then
  echo "number of \$kernel_dir/proc*_rho_kernel.bin files not equal to ${SEM_nproc_total}"
  echo "====== kernel process job FAILED: $event_id"
  exit 1
fi

echo ====== convert cijkl kernel to VTI kernel

mesh_dir=${SEM_iter_dir}/mesh/DATABASES_MPI
kernel_dir=${event_dir}/output_kernel/kernel
out_dir=${event_dir}/kernel/GLL
[ -d \${out_dir} ] && rm -rf \${out_dir}
mkdir -p \$out_dir

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_kernel_reparameterization.py \\
  --nproc ${SEM_nproc_total} \\
  --reference_dir ${SEM_reference_model_dir} \\
  --model_dir ${SEM_iter_dir}/model_initial \\
  --kernel_dir \${kernel_dir} \\
  --out_dir \${out_dir} \\
  --type ${SEM_parameterization_type} \\
  --mesh_dir \${mesh_dir}

# check if all files are created successfully
n=\$(ls -1 \$out_dir/*.bin | wc -l)
n0=\$(( ${SEM_nproc_total} * 6 )) # 6: number of kernel files for each process (alpha,beta,phi,xi,eta,rho)
if [ "\$n" -ne "\$n0" ]; then
  echo "====== kernel process job FAILED: $event_id"
  exit 1
fi

# remove forward saved frames adn adjoint source files to save disk space
chmod u+w $event_dir/forward_saved_frames
rm -rf $event_dir/forward_saved_frames
rm -rf $event_dir/SEM/*.adj
rm -rf \${kernel_dir} # remove kernel directory to save disk space

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo
EOF


#====== mask source/receiver 
cat <<EOF > $mask_job
#!/bin/bash
#SBATCH -J ${event_id}.mask
#SBATCH -o $mask_job.o%j
#SBATCH ${SLURM_args_kernel_process}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

echo ====== create source/receiver mask

mesh_dir=${SEM_iter_dir}/mesh
out_dir=${event_dir}/kernel/mask
[ -d \${out_dir} ] && rm -rf \${out_dir}
mkdir -p \$out_dir

echo "x,y,z,sigma_km" > \${out_dir}/mask.lst

awk 'NR==6{print \$1 "," \$2 "," \$3 "," a}' a="$SEM_kernel_mask_source_sigma_km" \\
  ${event_dir}/output_kernel/source.vtk \\
  >> \${out_dir}/mask.lst

awk 'NR>=6&&NF==3{print \$1 "," \$2 "," \$3 "," a}' a=$SEM_kernel_mask_receiver_sigma_km \\
  ${event_dir}/output_kernel/receiver.vtk \\
  >> \${out_dir}/mask.lst

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_make_gaussian_mask.py \\
  ${SEM_nproc_total} \\
  \${mesh_dir}/DATABASES_MPI \\
  \${out_dir}/mask.lst \\
  \${out_dir}

echo ====== mask kernels

out_dir=${event_dir}/kernel/GLL_mask
[ -d \${out_dir} ] && rm -rf \${out_dir}
mkdir -p \$out_dir

for model in ${SEM_model_names[@]}
do
  tag=\${model}${SEM_kernel_tag}
  ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_gll_math.py \\
    --nproc ${SEM_nproc_total} \\
    --model_dirs ${event_dir}/kernel/GLL  ${event_dir}/kernel/mask \\
    --model_tags \${tag} mask \\
    --math_expr "v[0]*v[1]" \\
    --out_dir \${out_dir} \\
    --out_tag \${tag} \\
    --overwrite_ok
done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo
EOF


# #====== kernel processing
# cat <<EOF > $process_job
# #!/bin/bash
# #SBATCH -J ${event_id}.process
# #SBATCH -o $process_job.o%j
# #SBATCH ${SLURM_args_kernel_process}
# 
# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
# echo
# 
# kernel_dir=${event_dir}/output_kernel/kernel
# 
# echo ====== check if simulation finished successfully
# 
# n=\$(ls -1 \$kernel_dir/proc*_cijkl_kernel.bin | wc -l)
# if [ "\$n" -ne "${SEM_nproc_total}" ]; then
#   echo "number of \$kernel_dir/proc*_cijkl_kernel.bin files not equal to ${SEM_nproc_total}"
#   echo "====== kernel process job FAILED: $event_id"
#   exit 1
# fi
# 
# n=\$(ls -1 \$kernel_dir/proc*_rho_kernel.bin | wc -l)
# if [ "\$n" -ne "${SEM_nproc_total}" ]; then
#   echo "number of \$kernel_dir/proc*_rho_kernel.bin files not equal to ${SEM_nproc_total}"
#   echo "====== kernel process job FAILED: $event_id"
#   exit 1
# fi
# 
# echo ====== convert cijkl kernel to VTI kernel
# 
# mesh_dir=${SEM_iter_dir}/mesh/DATABASES_MPI
# kernel_dir=${event_dir}/output_kernel/kernel
# out_dir=${event_dir}/kernel/GLL
# [ -d \${out_dir} ] && rm -rf \${out_dir}
# mkdir -p \$out_dir
# 
# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_kernel_reparameterization.py \\
#   --nproc ${SEM_nproc_total} \\
#   --reference_dir ${SEM_reference_model_dir} \\
#   --model_dir ${SEM_iter_dir}/model_initial \\
#   --kernel_dir \${kernel_dir} \\
#   --out_dir \${out_dir} \\
#   --type ${SEM_parameterization_type} \\
#   --mesh_dir \${mesh_dir}
# 
# # check if all files are created successfully
# n=\$(ls -1 \$out_dir/*.bin | wc -l)
# n0=\$(( ${SEM_nproc_total} * 6 )) # 6: number of kernel files for each process (alpha,beta,phi,xi,eta,rho)
# if [ "\$n" -ne "\$n0" ]; then
#   echo "====== kernel process job FAILED: $event_id"
#   exit 1
# fi
# 
# # remove forward saved frames adn adjoint source files to save disk space
# chmod u+w $event_dir/forward_saved_frames
# rm -rf $event_dir/forward_saved_frames
# rm -rf $event_dir/SEM/*.adj
# rm -rf \${kernel_dir} # remove kernel directory to save disk space
# 
# echo ====== create source/receiver mask
# 
# mesh_dir=${SEM_iter_dir}/mesh
# out_dir=${event_dir}/kernel/mask
# [ -d \${out_dir} ] && rm -rf \${out_dir}
# mkdir -p \$out_dir
# 
# echo "x,y,z,sigma_km" > \${out_dir}/mask.lst
# 
# awk 'NR==6{print \$1 "," \$2 "," \$3 "," a}' a="$SEM_kernel_mask_source_sigma_km" \\
#   ${event_dir}/output_kernel/source.vtk \\
#   >> \${out_dir}/mask.lst
# 
# awk 'NR>=6&&NF==3{print \$1 "," \$2 "," \$3 "," a}' a=$SEM_kernel_mask_receiver_sigma_km \\
#   ${event_dir}/output_kernel/receiver.vtk \\
#   >> \${out_dir}/mask.lst
# 
# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_make_gaussian_mask.py \\
#   ${SEM_nproc_total} \\
#   \${mesh_dir}/DATABASES_MPI \\
#   \${out_dir}/mask.lst \\
#   \${out_dir}
# 
# echo ====== mask kernels
# 
# out_dir=${event_dir}/kernel/GLL_mask
# [ -d \${out_dir} ] && rm -rf \${out_dir}
# mkdir -p \$out_dir
# 
# for model in ${SEM_model_names[@]}
# do
#   tag=\${model}${SEM_kernel_tag}
#   ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_gll_math.py \\
#     --nproc ${SEM_nproc_total} \\
#     --model_dirs ${event_dir}/kernel/GLL  ${event_dir}/kernel/mask \\
#     --model_tags \${tag} mask \\
#     --math_expr "v[0]*v[1]" \\
#     --out_dir \${out_dir} \\
#     --out_tag \${tag} \\
#     --overwrite_ok
# done
# 
# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
# echo
# EOF


# #====== kernel clipping   
# cat <<EOF > $threshold_job
# #!/bin/bash
# #SBATCH -J ${event_id}.threshold
# #SBATCH -o $threshold_job.o%j
# #SBATCH ${SLURM_args_kernel_threshold}
# 
# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
# echo
# 
# echo
# echo ====== convert cijkl kernel to VTI kernel
# echo
# 
# mesh_dir=${SEM_iter_dir}/mesh/DATABASES_MPI
# kernel_dir=${event_dir}/output_kernel/kernel
# out_dir=${event_dir}/kernel/GLL
# [ -d \${out_dir} ] && rm -rf \${out_dir}
# mkdir -p \$out_dir
# 
# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_kernel_reparameterization.py \\
#   --nproc ${SEM_nproc_total} \\
#   --reference_dir ${SEM_reference_model_dir} \\
#   --model_dir ${SEM_iter_dir}/model_initial \\
#   --kernel_dir \${kernel_dir} \\
#   --out_dir \${out_dir} \\
#   --type ${SEM_parameterization_type} \\
#   --mesh_dir \${mesh_dir}
# 
# # rm -rf \${kernel_dir} # remove kernel directory to save disk space
# 
# echo
# echo ====== apply threshold
# echo
# 
# kernel_dir=${event_dir}/kernel/GLL
# out_dir=${event_dir}/kernel/GLL_threshold
# [ -d \${out_dir} ] && rm -rf \${out_dir}
# mkdir -p \$out_dir
# 
# for model_name in ${SEM_model_names[@]}
# do
#   ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_gll_histogram.py  \\
#     ${SEM_nproc_total} \\
#     \${mesh_dir} \\
#     \${kernel_dir} \\
#     \${model_name}${SEM_kernel_tag} \\
#     --nbin ${SEM_histogram_nbin} \\
#     --exponential_base ${SEM_histogram_exponential_base} \\
#     --out_hist \${out_dir}/histogram_\${model_name}.txt \\
#     --cdf_threshold ${SEM_cdf_threshold} \\
#     --out_dir \${out_dir} 
# done
# 
# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
# echo
# EOF


#====== plot misfit and waveforms
cat <<EOF > $plot_job
#!/bin/bash
#SBATCH -J ${event_id}.plot
#SBATCH -o $plot_job.o%j
#SBATCH ${SLURM_args_plot}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

cd $event_dir

if [ -d "$figure_dir" ]
then
  chmod -R u+w $figure_dir
  rm -rf $figure_dir
fi
mkdir -p $figure_dir

$SEM_python_exec $SEM_utils_dir/misfit/plot.py \\
  $db_file $figure_dir \\
  --nproc=\${SLURM_NTASKS} \\
  --title="${SEM_stage_tag}_s${SEM_stage_num}i${SEM_iter_num}"

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF


#====== perturb: forward simulation of perturbed model
cat <<EOF > $perturb_job
#!/bin/bash
#SBATCH -J ${event_id}.perturb
#SBATCH -o $perturb_job.o%j
#SBATCH ${SLURM_args_forward}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

cd $event_dir/DATA
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file

#for tag in dvp dvs
for dm_tag in ${SEM_perturb_group_names[@]}
do

  mesh_perturb_dir=$SEM_iter_dir/mesh_perturb_\${dm_tag}
  out_dir=output_perturb_\${dm_tag}

  rm -rf $event_dir/DATABASES_MPI
  mkdir $event_dir/DATABASES_MPI
  ln -s \$mesh_perturb_dir/DATABASES_MPI/*.bin $event_dir/DATABASES_MPI

  cd $event_dir

  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES

  cp \$mesh_perturb_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES

  ${SLURM_mpiexec} $SEM_build_dir/bin/xspecfem3D

  # check if simulation finished successfully
  if grep -qF "End of the simulation" \$out_dir/output_solver.txt; then
    echo "====== forward perturb job SUCCESS: $event_id \${dm_tag}"
  else
    echo "====== forward perturb job FAILED: $event_id \${dm_tag}"
    exit 1
  fi

  mkdir $event_dir/\$out_dir/sac
  mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo
EOF


#====== search: cc linearized seismograms for chosen step sizes
cat <<EOF > $search_job
#!/bin/bash
#SBATCH -J ${event_id}.search
#SBATCH -o $search_job.o%j
#SBATCH ${SLURM_args_search}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

echo
echo "Read synthetic waveforms of perturbed model [\$(date -Is)]"
echo

# for tag in dxs dmt
for dm_tag in ${SEM_perturb_group_names[@]}
do
  $SEM_python_exec $SEM_utils_dir/misfit/read_perturbed_syn.py \\
    $db_file \\
    $event_dir/output_perturb_\${dm_tag}/sac \\
    \${dm_tag}
done

echo
echo "grid search [\$(date -Is)]"
echo

export OPENBLAS_NUM_THREADS=1

$SEM_python_exec $SEM_utils_dir/misfit/grid_search_structure.py \\
  $db_file \\
  --dm_tags ${SEM_perturb_group_names[@]} \\
  --dm_steps ${SEM_search_step_sizes[@]} \\
  --nproc=\$SLURM_NPROCS

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

exit -1

#///////////////////////// Hessian simulation

#====== hess_syn: forward simulation of randomly perturbed model
cat <<EOF > $hess_syn_job
#!/bin/bash
#SBATCH -J ${event_id}.hess_syn
#SBATCH -o $hess_syn_job.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_hess_forward
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir/DATA
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .true./" Par_file

for tag in ${hess_model_names}
do

  echo ====== \$tag

  out_dir=output_syn_\${tag}

  rm -rf $event_dir/DATABASES_MPI
  mkdir $event_dir/DATABASES_MPI
  ln -s $iter_dir/mesh_\${tag}/DATABASES_MPI/*.bin $event_dir/DATABASES_MPI

  cd $event_dir

  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES

  cp $iter_dir/mesh_\${tag}/OUTPUT_FILES/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES

  ${slurm_mpiexec} $sem_build_dir/bin/xspecfem3D

  mkdir $event_dir/\$out_dir/sac
  mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== hess_misfit
cat <<EOF > $hess_misfit_job
#!/bin/bash
#SBATCH -J ${event_id}.hess_misfit
#SBATCH -o $hess_misfit_job.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_hess_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir

for tag in ${hess_model_names}
do

  echo ====== \$tag

  # make adjoint source
  rm -rf $event_dir/adj_kernel_\${tag}
  mkdir -p $event_dir/adj_kernel_\${tag}

  $utils_dir/output_adj_for_perturbed_waveform.py \
    $misfit_par_file \
    $db_file \
    $event_dir/output_syn_\${tag}/sac \
    $event_dir/adj_kernel_\${tag}

  # make STATIONS_ADJOINT
  cd $event_dir/adj_kernel_\${tag}
  ls *Z.adj | sed 's/..Z\.adj$//' |\
    awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
  grep -f $event_dir/adj_kernel_\${tag}/grep_pattern $event_dir/DATA/STATIONS \
    > $event_dir/adj_kernel_\${tag}/STATIONS_ADJOINT

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== kernel simulation for randomly perturbed model
cat <<EOF > $hess_kernel_job
#!/bin/bash
#SBATCH -J ${event_id}.hess_kernel
#SBATCH -o $hess_kernel_job.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_hess_adjoint
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

for tag in ${hess_model_names}
do

  echo ====== \${tag}

  out_dir=output_kernel_\${tag}

  cd $event_dir/DATA
  chmod u+w Par_file
  sed -i "/^SIMULATION_TYPE/s/=.*/= 3/" Par_file
  sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
  sed -i "/^ANISOTROPIC_KL/s/=.*/= .true./" Par_file
  sed -i "/^SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
  sed -i "/^APPROXIMATE_HESS_KL/s/=.*/= .false./" Par_file

  cp -f $event_dir/adj_kernel_\${tag}/STATIONS_ADJOINT $event_dir/DATA/
  rm -rf $event_dir/SEM
  ln -s $event_dir/adj_kernel_\${tag} $event_dir/SEM

  cd $event_dir

  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES

  cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES

  cd $event_dir
  ${slurm_mpiexec} $sem_build_dir/bin/xspecfem3D

  mkdir $event_dir/\$out_dir/sac
  mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

  mkdir $event_dir/\$out_dir/kernel
  mv $event_dir/DATABASES_MPI/*reg1_cijkl_kernel.bin $event_dir/\$out_dir/kernel/
  mv $event_dir/DATABASES_MPI/*reg1_rho_kernel.bin $event_dir/\$out_dir/kernel/

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF
