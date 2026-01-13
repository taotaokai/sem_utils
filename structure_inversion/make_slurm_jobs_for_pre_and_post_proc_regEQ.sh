#!/bin/bash
# Make jobs files for slurm
# structure inversion with regional earthquake

control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}

#====== source control_file
source $control_file
event_list=$(readlink -f $event_list)

#====== define directories
#mesh_dir=$SEM_iter_dir/mesh # DATABASES_MPI/proc*_reg1_solver_data.bin
#hess_dir=$iter_dir/hess_sum
#mesh_perturb_dir=$iter_dir/mesh_perturb # DATABASES_MPI/proc*_reg1_solver_data.bin
#model_dir=$iter_dir/model # proc*_reg1_vph,vpv,vsv,vsh,eta,rho,dlnvs,kappa,eta,gamma.bin
#model_random_dir=$iter_dir/model_random # for model with small random perturbations
#kernel_dir=$iter_dir/kernel # DATABASES_MPI/proc*_reg1_solver_data.bin
slurm_dir=$SEM_iter_dir/slurm
mkdir -p $slurm_dir
#utils_dir=$sem_utils_dir/utils/misfit_v0 # sem_utils/utils/misfit_v0

#====== define job names
mesh_job=$slurm_dir/mesh.job

kernel_sum_job=$slurm_dir/kernel_sum.job
kernel_precond_job=$slurm_dir/kernel_precond.job

dmodel_job=$slurm_dir/dmodel.job
model_perturb_job=$slurm_dir/model_perturb.job
mesh_perturb_job=$slurm_dir/mesh_perturb.job
model_update_job=$slurm_dir/model_update.job

model_perturb_random_job=$slurm_dir/model_perturb_random.job

#====== mesh
cat <<EOF > $mesh_job
#!/bin/bash
#SBATCH -J mesh
#SBATCH -o ${mesh_job}.o%j
#SBATCH ${SLURM_args_mesh}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

mesh_dir=${SEM_iter_dir}/mesh
model_dir=${SEM_iter_dir}/model_initial

if [ ! -d "\$model_dir" ]
then
  echo "[ERROR] \$model_dir does not exist!"
  exit -1
fi

#rm -rf \$mesh_dir
mkdir -p \$mesh_dir

cd \$mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

cd \$mesh_dir/DATA
chmod u+w *
ln -sf $SEM_build_dir/DATA/* \$mesh_dir/DATA/
rm Par_file GLL CMTSOLUTION
ln -sf \$model_dir GLL
cp -L $SEM_config_dir/DATA/Par_file .
cp -L $SEM_config_dir/DATA/CMTSOLUTION .
cp -L Par_file \$mesh_dir/OUTPUT_FILES/

sed -i '/^MODEL/s/=[[:space:]]*[^[:space:]_#]*/= GLL/' \$mesh_dir/DATA/Par_file
sed -i '/^SAVE_MESH_FILES/s/=.*/= .false./' \$mesh_dir/DATA/Par_file
#sed -i "/^USE_ECEF_CMTSOLUTION/s/=.*/= .false./" \$mesh_dir/DATA/Par_file
#sed -i "/^USE_FORCE_POINT_SOURCE/s/=.*/= .false./" \$mesh_dir/DATA/Par_file

cd \$mesh_dir
${SLURM_mpiexec} $SEM_build_dir/bin/xmeshfem3D

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

#====== kernel_sum
cat <<EOF > $kernel_sum_job
#!/bin/bash
#SBATCH -J kernel_sum
#SBATCH -o ${kernel_sum_job}.o%j
#SBATCH ${SLURM_args_kernel_sum}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

echo "====== sum up all event kernels"

kernel_dir=${SEM_iter_dir}/kernel_sum
mkdir -p \$kernel_dir

awk 'NF&&\$1!~/#/{printf "%s/events/%s/kernel/GLL_threshold\\n", a,\$1}' \\
  a="$SEM_iter_dir" $event_list > \\
  \${kernel_dir}/event_kernel_dir.list

for model_name in ${SEM_model_names[@]}
do
  ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_sum.py \\
    ${SEM_nproc_total} \\
    \${kernel_dir}/event_kernel_dir.list \\
    \${model_name}${SEM_kernel_tag} \\
    \${kernel_dir} \\
    --ncomp=1
done

# awk 'NF&&\$1!~/#/{printf "%s/events/%s/output_kernel/kernel\\n", a,\$1}' \\
#   a="$SEM_iter_dir" $event_list > \\
#   \${kernel_dir}/event_kernel_dir.list
# 
# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_sum.py \\
#   ${SEM_nproc_total} \\
#   \${kernel_dir}/event_kernel_dir.list \\
#   cijkl_kernel \\
#   \${kernel_dir} \\
#   --mask_tag=mask --ncomp=21
# 
# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_sum.py \\
#   ${SEM_nproc_total} \\
#   \${kernel_dir}/event_kernel_dir.list \\
#   rho_kernel \\
#   \${kernel_dir} \\
#   --mask_tag=mask --ncomp=1

# delete event kernel directories to save disk space
# for event_kernel_dir in \$(cat \${kernel_dir}/event_kernel_dir.list)
# do
#   rm -rf \$event_kernel_dir
# done

# echo "====== readuce cijkl,rho_kernel to VTI model parameters"
# 
# kernel_dir=${SEM_iter_dir}/kernel_sum
# 
# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_kernel_in_alpha_beta_phi_xi_eta_from_cijkl_rho.py \\
#   ${SEM_nproc_total} \\
#   ${SEM_iter_dir}/model_initial \\
#   ${SEM_reference_model_dir} \\
#   \${kernel_dir} \\
#   \${kernel_dir}

echo
echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

#====== kernel_precond
cat <<EOF > $kernel_precond_job
#!/bin/bash
#SBATCH -J kernel_precond
#SBATCH -o ${kernel_precond_job}.o%A_%a
#SBATCH --array 0-$(( ${#SEM_model_names[@]} - 1 ))
#SBATCH ${SLURM_args_kernel_precond}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

model_names=(${SEM_model_names[@]})

kernel_name=\${model_names[\${SLURM_ARRAY_TASK_ID}]}${SEM_kernel_tag}

echo "====== smooth kernel"

mesh_dir=${SEM_iter_dir}/mesh
kernel_dir=${SEM_iter_dir}/kernel_sum

out_dir=${SEM_iter_dir}/kernel_precond
mkdir -p \$out_dir

# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_smooth_diffusion_iso.py \\
#   ${SEM_nproc_total} \\
#   \${mesh_dir}/DATABASES_MPI \\
#   \${kernel_dir} \\
#   \${kernel_name} \\
#   ${SEM_kernel_smooth_iso_FWHM_km} \\
#   ${SEM_kernel_smooth_diffusion_niter} \\
#   \${out_dir} \\
#   --max_tolerance=1e-5

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_smooth_diffusion_VTI.py \\
  ${SEM_nproc_total} \\
  \${mesh_dir}/DATABASES_MPI \\
  \${kernel_dir} \\
  \${kernel_name} \\
  \${out_dir} \\
  --horizontal_length ${SEM_kernel_smooth_horizontal_length_km} \\
  --vertical_length ${SEM_kernel_smooth_vertical_length_km} \\
  --nstep ${SEM_kernel_smooth_diffusion_niter} \\
  --max_tolerance=1e-3

echo
echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

#TODO scale kernel_precond ?


#====== search direction
cat <<EOF > $dmodel_job
#!/bin/bash
#SBATCH -J dmodel
#SBATCH -o ${dmodel_job}.o%j
#SBATCH ${SLURM_args_dmodel}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

out_dir=${SEM_iter_dir}/dmodel
mkdir -p \${out_dir}

# get search direction

if [ "$SEM_iter_num" -eq 0 ]
then

  # first iteration use scaled gradient as search direction
  ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_scale2.py \\
    ${SEM_nproc_total} \\
    ${SEM_iter_dir}/kernel_precond \\
    \${out_dir} \\
    --models ${SEM_model_names[@]} \\
    --in_tag ${SEM_kernel_tag} \\
    --out_tag ${SEM_dmodel_tag} \\
    --scaled_amplitude=${SEM_dmodel_scale}

else

  # search direction by CG method: Hestenes and Stiefel (1952)
  ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_cg_maximize2.py \\
    ${SEM_nproc_total} \\
    ${SEM_iter_dir}/mesh/DATABASES_MPI \\
    ${SEM_prev_iter_dir}/model_initial \\
    ${SEM_iter_dir}/model_initial \\
    ${SEM_prev_iter_dir}/kernel_sum \\
    ${SEM_iter_dir}/kernel_sum \\
    ${SEM_iter_dir}/kernel_precond \\
    \${out_dir} \\
    --models ${SEM_model_names[@]} \\
    --kernel_tag ${SEM_kernel_tag} \\
    --dmodel_tag ${SEM_dmodel_tag} \\
    --scaled_amplitude=${SEM_dmodel_scale} \\
    --cg_type=${SEM_cg_type}

fi

echo
echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF


#====== model_perturb
cat <<EOF > $model_perturb_job
#!/bin/bash
#SBATCH -J model_perturb
#SBATCH -o ${model_perturb_job}.o%j
#SBATCH ${SLURM_args_model_perturb}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

#====== model perturbation for each group

group_names=(${SEM_perturb_group_names[@]})
model_groups=(${SEM_perturb_model_groups[@]})

ngroup=\${#group_names[@]}
for ((i=0; i<\$ngroup; i++))
do

  dm_tag=\${group_names[\$i]} 
  out_dir=$SEM_iter_dir/model_perturb_\${dm_tag}
  [ -e \$out_dir ] && rm -rf \$out_dir
  mkdir -p \$out_dir

  ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_perturb_groups.py \\
    ${SEM_nproc_total} \\
    $SEM_iter_dir/model_initial \\
    $SEM_iter_dir/dmodel \\
    \$out_dir \\
    --model_groups \${model_groups[\$i]} \\
    --group_scales 1.0 \\
    --dmodel_tag ${SEM_dmodel_tag} \\
    --method "absolute" # \\
    # --models ${SEM_model_names[@]} \\
    # --model_min_values ${SEM_model_min_values[@]} \\
    # --model_max_values ${SEM_model_max_values[@]}

done

#   ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_perturb.py \\
#     ${SEM_nproc_total} \\
#     $SEM_iter_dir/model_initial \\
#     $SEM_iter_dir/dmodel \\
#     \$out_dir \\
#     --model_tags \${perturb_model_tag_groups[\$i]//,/ } \\
#     --dmodel_tags \${perturb_dmodel_tag_groups[\$i]//,/ } \\
#     --min_values \${perturb_min_values[\$i]//,/ } \\
#     --max_values \${perturb_max_values[\$i]//,/ } \\
#     --scale 1.0 \\
#     --method "absolute"
# 
# done

#====== convert model from alpha,beta,xi,phi to vpv,vph,vsv,vsh

# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh.py \\
#   ${SEM_nproc_total} \\
#   ${SEM_reference_model_dir} \\
#   \$out_dir \\
#   \$out_dir

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_model_reparameterization.py \\
  ${SEM_nproc_total} \\
  ${SEM_reference_model_dir} \\
  \$out_dir \\
  \$out_dir \\
  --reverse \\
  --type ${SEM_parameterization_type}

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF


#====== mesh_perturb
cat <<EOF > $mesh_perturb_job
#!/bin/bash
#SBATCH -J mesh_perturb
#SBATCH -o ${mesh_perturb_job}.o%j
#SBATCH ${SLURM_args_mesh_perturb}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

group_names=(${SEM_perturb_group_names[@]})
ngroup=${#SEM_perturb_group_names[@]}

for ((i=0; i<\$ngroup; i++))
do

  dm_tag=\${group_names[\$i]} 

  model_dir=${SEM_iter_dir}/model_perturb_\${dm_tag}
  if [ ! -d "\$model_dir" ]
  then
    echo "[ERROR] \$model_dir does not exist!"
    exit -1
  fi

  mesh_dir=${SEM_iter_dir}/mesh_perturb_\${dm_tag}
  [ -e \$mesh_dir ] && rm -rf \$mesh_dir
  mkdir -p \$mesh_dir

  cd \$mesh_dir
  mkdir DATA DATABASES_MPI OUTPUT_FILES

  cd \$mesh_dir/DATA
  chmod u+w *
  ln -sf $SEM_build_dir/DATA/* \$mesh_dir/DATA/
  rm Par_file GLL CMTSOLUTION

  mkdir GLL
  ln -s -t GLL ${SEM_iter_dir}/model_initial/*.bin # first link initial model files
  ln -sf -t GLL \$model_dir/*.bin  # -f force override with links to new model files 

  cp -L $SEM_config_dir/DATA/Par_file .
  cp -L $SEM_config_dir/DATA/CMTSOLUTION .
  cp -L Par_file \$mesh_dir/OUTPUT_FILES/

  sed -i '/^MODEL/s/=[[:space:]]*[^[:space:]_#]*/= GLL/' \$mesh_dir/DATA/Par_file
  sed -i '/^SAVE_MESH_FILES/s/=.*/= .false./' \$mesh_dir/DATA/Par_file

  cd \$mesh_dir
  ${SLURM_mpiexec} $SEM_build_dir/bin/xmeshfem3D

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF


#====== model_update
cat <<EOF > $model_update_job
#!/bin/bash
#SBATCH -J model_update
#SBATCH -o ${model_update_job}.o%j
#SBATCH ${SLURM_args_model_update}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

echo ====== collect grid search results from all events

out_dir=${SEM_iter_dir}/line_search
mkdir -p \${out_dir}

awk 'NF&&\$1!~/#/{printf "%s/events/%s/misfit/misfit.h5\\n", a,\$1}' \\
  a="$SEM_iter_dir" $event_list > \\
  \${out_dir}/misfit_h5file.list

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/structure_inversion/grid_search.py \\
  \${out_dir}/misfit_h5file.list \\
  --out_nc \${out_dir}/grid_search.nc \\
  --out_figure \${out_dir}/grid_search.pdf \\
  --out_txt \${out_dir}/grid_search.txt

echo ====== apply model update

if [ -d "$SEM_iter_dir/model_updated" ]
then
  rm -rf $SEM_iter_dir/model_updated
fi
mkdir -p $SEM_iter_dir/model_updated

# copy initial model in case some model parameters are not inverted, e.g., qmu, rho 
cp -L $SEM_iter_dir/model_initial/*.bin $SEM_iter_dir/model_updated/
chmod -R u+w $SEM_iter_dir/model_updated

group_names=(${SEM_perturb_group_names[@]})
model_groups=(${SEM_perturb_model_groups[@]})
ngroup=\${#group_names[@]}

# get optimal step length for each model search direction
group_scales=()
for ((i=0; i<\$ngroup; i++))
do
  dm_tag=\${group_names[\$i]} 
  opt_dm_scale=\$(grep \$dm_tag \${out_dir}/grid_search.txt | awk '{print \$NF}')
  group_scales[\$i]=\$opt_dm_scale
  echo "# apply model update for [\${dm_tag}] with step length [\${opt_dm_scale}]"
done

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_perturb_groups.py \\
  ${SEM_nproc_total} \\
  $SEM_iter_dir/model_initial \\
  $SEM_iter_dir/dmodel \\
  $SEM_iter_dir/model_updated \\
  --model_groups \${model_groups[@]} \\
  --group_scales \${group_scales[@]} \\
  --dmodel_tag ${SEM_dmodel_tag} \\
  --method "absolute" \\
  --models ${SEM_model_names[@]} \\
  --model_min_values ${SEM_model_min_values[@]} \\
  --model_max_values ${SEM_model_max_values[@]}


# ngroup=\${#perturb_group_names[@]}
# 
# for ((i=0; i<\$ngroup; i++))
# do
# 
#   dm_tag=\${perturb_group_names[\$i]} 
# 
#   opt_dm_scale=\$(grep \$dm_tag \${out_dir}/grid_search.txt | awk '{print \$NF}')
# 
#   echo "# apply model update for [\${dm_tag}] with step length [\${opt_dm_scale}]"
# 
#   ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_perturb.py \\
#     ${SEM_nproc_total} \\
#     $SEM_iter_dir/model_initial \\
#     $SEM_iter_dir/dmodel \\
#     $SEM_iter_dir/model_updated \\
#     --model_tags \${perturb_model_tag_groups[\$i]//,/ } \\
#     --dmodel_tags \${perturb_dmodel_tag_groups[\$i]//,/ } \\
#     --min_values \${perturb_min_values[\$i]//,/ } \\
#     --max_values \${perturb_max_values[\$i]//,/ } \\
#     --scale \$opt_dm_scale \\
#     --method "absolute"
# 
# done

# ${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh.py \\
#   ${SEM_nproc_total} \\
#   ${SEM_reference_model_dir} \\
#   $SEM_iter_dir/model_updated \\
#   $SEM_iter_dir/model_updated \\
#   $model_update_min_alpha $model_update_max_alpha \\
#   $model_update_min_beta  $model_update_max_beta  \\
#   $model_update_min_phi   $model_update_max_phi   \\
#   $model_update_min_xi    $model_update_max_xi    \\
#   $model_update_min_eta   $model_update_max_eta

echo ====== convert back to vph,vpv,vsv,vsh

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_model_reparameterization.py \\
  ${SEM_nproc_total} \\
  ${SEM_reference_model_dir} \\
  $SEM_iter_dir/model_updated \\
  $SEM_iter_dir/model_updated \\
  --reverse \\
  --type ${SEM_parameterization_type}

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF

#====== random model perturbation
cat <<EOF > $model_perturb_random_job
#!/bin/bash
#SBATCH -J model_perturb_random
#SBATCH -o ${model_perturb_random_job}.o%j
#SBATCH ${SLURM_args_model_perturb}

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_gll_random_perturb.py \\
  ${SEM_nproc_total} \\
  ${SEM_iter_dir}/model_initial \\
  ${SEM_iter_dir}/model_perturb_random \\
  --model_tags ${SEM_model_names[@]} \\
  --amplitudes ${SEM_hessian_perturb_amplitude} \\
  --method "absolute"

# convert back to vpv,vph,vsv,vsh
${SLURM_mpiexec} ${SEM_python_exec} $SEM_utils_dir/meshfem3d/sem_VTI_model_reparameterization.py \\
  ${SEM_nproc_total} \\
  ${SEM_reference_model_dir} \\
  $SEM_iter_dir/model_perturb_random \\
  $SEM_iter_dir/model_perturb_random \\
  --reverse \\
  --type ${SEM_parameterization_type}

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date -Is)]"
echo

EOF



exit -1

# #====== kernel_convert_mask
# cat <<EOF > $kernel_convert_mask_job
# #!/bin/bash
# #SBATCH -J kernel_convert_mask
# #SBATCH -o ${kernel_convert_mask_job}.o%j
# #SBATCH ${slurm_args_kernel_convert_mask}
# 
# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -I)]"
# echo
# 
# echo "====== convert cijkl,rho_kernel to tiso kernel with masks"
# 
# for event_id in \$(awk 'NF&&\$1!~/#/{print \$1}' $event_list)
# do
#   echo "------ \$event_id"
#   
#   event_dir=${iter_dir}/events/\$event_id
#   ker_dir=\${event_dir}/output_kernel/kernel  # input files: cijkl_kernel, rho_kernel
#   out_dir=\${event_dir}/output_kernel/tiso_kernel_masked
#   mkdir -p \$out_dir
# 
#   awk 'NR==6{print \$0, a}' a="$sem_kernel_mask_source_sigma_km" \\
#     \${event_dir}/output_kernel/source.vtk \\
#     > \${out_dir}/mask.lst
# 
#   awk 'NR>=6&&NF==3{print \$0, a}' a=$sem_kernel_mask_receiver_sigma_km \\
#     \${event_dir}/output_kernel/receiver.vtk \\
#     >> \${out_dir}/mask.lst
# 
#   ${slurm_mpiexec} ${python_exec} $sem_utils_dir/meshfem3d/_from_cijkl_rho.py \\
#     ${sem_nproc_total} \\
#     \${event_dir}/DATABASES_MPI \\
#     \${ker_dir} \\
#     \${out_dir} \\
#     --with_mask \\
#     --mesh_dir \${event_dir}/DATABASES_MPI \\
#     --mask_list \${out_dir}/mask.lst
# done
# 
# echo
# echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date -I)]"
# echo
# 
# EOF

# #====== kernel_sum_smooth_precond
# cat <<EOF > $kernel_sum_smooth_precond_job
# #!/bin/bash
# #SBATCH -J kernel_sum_smooth_precond
# #SBATCH -o ${kernel_sum_smooth_precond_job}.o%A_%a
# #SBATCH --array 0-$(( ${#sem_kernel_tags[@]} - 1 ))
# #SBATCH ${slurm_args_kernel_sum_smooth_precond}
# 
# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date -I)]"
# echo
# 
# sem_kernel_tags=(${sem_kernel_tags[@]})
# 
# kernel_tag=\${sem_kernel_tags[\${SLURM_ARRAY_TASK_ID}]}
# 
# echo "====== sum up all event kernels"
# 
# kernel_dir=${iter_dir}/kernel
# mkdir -p \$kernel_dir
# 
# awk 'NF&&\$1!~/#/{printf "%s/events/%s/output_kernel/kernel\\n", a,\$1}' \\
#   a="$iter_dir" $event_list > \\
#   \${kernel_dir}/event_kernel_dir.list
# 
# ${slurm_mpiexec} ${python_exec} $sem_utils_dir/meshfem3d/sem_sum.py \\
#   ${sem_nproc_total} \\
#   \${kernel_dir}/event_kernel_dir.list \\
#   "\${kernel_tag}" \\
#   "\${kernel_dir}"
# 
# echo "====== smooth kernel"
# 
# out_dir=${iter_dir}/kernel_smooth_precond
# mkdir -p \$out_dir
# 
# ${slurm_mpiexec} ${python_exec} $sem_utils_dir/meshfem3d/sem_smooth_diffusion_iso.py \\
#   ${sem_nproc_total} \\
#   ${mesh_dir} \\
#   \${kernel_dir} \\
#   \${kernel_tag} \\
#   ${sem_kernel_smooth_iso_FWHM_km} \\
#   ${sem_kernel_smooth_diffusion_niter} \\
#   \${out_dir} \\
#   --max_tolerance=1e-5
# 
# echo
# echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date -I)]"
# echo
# 
# EOF


# #====== kernel_sum
# cat <<EOF > $kernel_sum_job
# #!/bin/bash
# #SBATCH -J kernel_sum
# #SBATCH -o ${kernel_sum_job}.o%j
# #SBATCH ${slurm_args_mesh}

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# echo "====== make source and receiver mask [\$(date)]"

# for event_id in \$(awk -F"|" 'NF&&\$1!~/#/{print \$9}' $event_list)
# do
#   echo "------ \$event_id"
#   event_dir=$iter_dir/\$event_id
#   out_dir=\$event_dir/source_receiver_mask
#   mkdir \$out_dir

#   awk 'NR==6{print \$0, a}' a=$source_mask_1sigma_km \$event_dir/output_kernel/source.vtk > \$out_dir/source.xyz
#   awk 'NR>=6&&NF==3{print \$0, a}' a=$receiver_mask_1sigma_km \$event_dir/output_kernel/receiver.vtk >> \$out_dir/source.xyz
#   ${slurm_mpiexec} $sem_utils_dir/meshfem3D/sem_make_gaussian_mask \
#     $sem_nproc $mesh_dir/DATABASES_MPI \
#     \$out_dir/source.xyz \
#     \$out_dir

#   ln -s \$out_dir/*_mask.bin \$event_dir/output_kernel/kernel
# done

# echo "====== sum up event kernels [\$(date)]"

# awk -F"|" 'NF&&\$1!~/#/{printf "%s/%s/output_kernel/kernel\\n", a,\$9}' \
#   a="$iter_dir" $event_list > kernel_dir.list

# out_dir=$iter_dir/kernel
# mkdir \$out_dir

# echo ------ sum up cijkl,rho_kernel with source and receiver mask

# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_sum_event_kernels_cijkl \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   kernel_dir.list cijkl_kernel \
#   1 "mask" \
#   \$out_dir cijkl_kernel

# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_sum_event_kernels_1 \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   kernel_dir.list rho_kernel \
#   1 "mask" \
#   \$out_dir rho_kernel

# echo "------ convert cijkl,rho_kernel to aijkl,rhoprime kernel [\$(date)]"

# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \
#   $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \
#   \$out_dir \
#   \$out_dir

# echo "------ reduce aijkl_kernel to alpha,beta,phi,xi,eta_kernel [\$(date)]"

# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_aijkl_to_tiso_in_alpha_beta_phi_xi_eta \
#   $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \
#   \$out_dir \
#   \$out_dir

# echo ====== precondition kernel

# for kernel_tag in alpha beta phi xi eta
# do
#   ${slurm_mpiexec} $sem_utils_dir/bin/xsem_math \
#     $sem_nproc $mesh_dir/DATABASES_MPI \
#     \$out_dir \${kernel_tag}_kernel \
#     $precond_dir inv_hess_diag \
#     "mul" \
#     \$out_dir \${kernel_tag}_kernel_precond
# done

# model_tags=alpha_kernel_precond,beta_kernel_precond,phi_kernel_precond,xi_kernel_precond,eta_kernel_precond

# $slurm_mpiexec $sem_utils_dir/bin/xsem_smooth \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   \$out_dir \${model_tags} \
#   $kernel_smooth_1sigma_h $kernel_smooth_1sigma_v \
#   \$out_dir "_smooth"

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== pcg_dmodel
# cat <<EOF > $pcg_dmodel_job
# #!/bin/bash
# #SBATCH -J pcg_dmodel
# #SBATCH -o ${pcg_dmodel_job}.o%j
# #SBATCH -N 1
# #SBATCH -n $nproc_per_node
# #SBATCH -p $slurm_partition
# #SBATCH -t $slurm_timelimit_misfit
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# #------ preconditioned conjugate gradient
# if [ x$iter_num != x00 ]
# then

#   current_precond_kernel_names=alpha_kernel_precond_smooth,beta_kernel_precond_smooth,phi_kernel_precond_smooth,xi_kernel_precond_smooth,eta_kernel_precond_smooth
#   current_kernel_names=alpha_kernel,beta_kernel,phi_kernel,xi_kernel,eta_kernel
#   previous_kernel_names=alpha_kernel,beta_kernel,phi_kernel,xi_kernel,eta_kernel
#   previous_dmodel_names=alpha_dmodel,beta_dmodel,phi_dmodel,xi_dmodel,eta_dmodel
#   out_dmodel_names=alpha_dmodel,beta_dmodel,phi_dmodel,xi_dmodel,eta_dmodel

#   previous_kernel_dir=$prev_iter_dir/kernel
#   previous_dmodel_dir=$prev_iter_dir/model_updated

#   $slurm_mpiexec $sem_utils_dir/bin/xsem_pcg_dmodel_n \
#     $sem_nproc $mesh_dir/DATABASES_MPI \
#     $kernel_dir  \$current_precond_kernel_names \
#     $kernel_dir  \$current_kernel_names \
#     \$previous_kernel_dir \$previous_kernel_names \
#     \$previous_dmodel_dir \$previous_dmodel_names \
#     $cg_type \
#     $kernel_dir \$out_dmodel_names

# else

#   cd $kernel_dir
#   ls *_precond_smooth.bin | awk -F"_" '{printf "ln -s %s %s_%s_%s_dmodel.bin\\n",\$0,\$1,\$2,\$3}' > ln_sh
#   bash ln_sh

# fi

# #------ scale dmodel
# kernel_suffix="_dmodel"
# out_suffix="_dmodel_scale"

# $slurm_mpiexec $sem_utils_dir/bin/xsem_scale_dmodel_alpha_beta_phi_xi_eta \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   $kernel_dir \$kernel_suffix \
#   $dmodel_max_dlnvs \
#   $kernel_dir \$out_suffix

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== model_perturb
# cat <<EOF > $model_perturb_job
# #!/bin/bash
# #SBATCH -J model_perturb
# #SBATCH -o ${model_perturb_job}.o%j
# #SBATCH -N 1
# #SBATCH -n $nproc_per_node
# #SBATCH -p $slurm_partition
# #SBATCH -t $slurm_timelimit_misfit
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# out_dir=$iter_dir/model_perturb
# mkdir \$out_dir

# #====== add dmodel
# model_names=alpha,beta,phi,xi,eta
# dmodel_names=alpha_dmodel_scale,beta_dmodel_scale,phi_dmodel_scale,xi_dmodel_scale,eta_dmodel_scale
# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_math \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   $model_dir \$model_names \
#   $kernel_dir \$dmodel_names \
#   "add" \
#   \$out_dir \$model_names

# #====== convert to gll model

# # link reference model
# ln -sf $model_REF_dir/proc*_reg1_v[ps]0.bin \$out_dir
# ln -sf $model_REF_dir/proc*_reg1_rho0.bin \$out_dir

# $slurm_mpiexec $sem_utils_dir/bin/xsem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh \
#   $sem_nproc $mesh_dir/DATABASES_MPI \$out_dir  \$out_dir

# #====== scale rho to beta
# $slurm_mpiexec $sem_utils_dir/bin/xsem_gll_scale_rho_to_beta \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   \$out_dir \$out_dir $model_scale_rho_to_beta

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== mesh_perturb
# cat <<EOF > $mesh_perturb_job
# #!/bin/bash
# #SBATCH -J mesh_perturb
# #SBATCH -o ${mesh_perturb_job}.o%j
# #SBATCH -N $slurm_nnode
# #SBATCH -n $slurm_nproc
# #SBATCH -p $slurm_partition
# #SBATCH -t $slurm_timelimit_misfit
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# mesh_dir=$iter_dir/mesh_perturb
# model_dir=$iter_dir/model_perturb

# if [ ! -d "\$model_dir" ]
# then
#   echo "[ERROR] \$model_dir does not exist!"
#   exit -1
# fi

# #rm -rf \$mesh_dir
# mkdir -p \$mesh_dir

# cd \$mesh_dir
# mkdir DATA DATABASES_MPI OUTPUT_FILES

# cd \$mesh_dir/DATA
# ln -sf $sem_build_dir/DATA/* \$mesh_dir/DATA/
# rm Par_file GLL CMTSOLUTION
# ln -sf \$model_dir GLL
# cp -L $sem_config_dir/DATA/Par_file .
# cp -L $sem_config_dir/DATA/CMTSOLUTION .
# cp -L Par_file CMTSOLUTION \$mesh_dir/OUTPUT_FILES/

# cd \$mesh_dir
# sed -i "/^MODEL/s/=.*/= GLL/" \$mesh_dir/DATA/Par_file

# ${slurm_mpiexec} $sem_build_dir/bin/xmeshfem3D

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== model_update
# cat <<EOF > $model_update_job
# #!/bin/bash
# #SBATCH -J model_update
# #SBATCH -o ${model_update_job}.o%j
# #SBATCH -N 1
# #SBATCH -n $nproc_per_node
# #SBATCH -p $slurm_partition
# #SBATCH -t $slurm_timelimit_misfit
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# cd $iter_dir
# mkdir grid_search_dmodel
# ls */misfit/grid_search_dmodel.txt | awk -F"/" '{printf "cp %s grid_search_dmodel/%s.txt\\n",\$0,\$1}' | bash

# cd $iter_dir/grid_search_dmodel
# ls *.txt > list
# $sem_utils_dir/utils/misfit_v0/plot_grid_search_1d.py list > grid_search.out
# step_length=\$(grep "optimal step_length" $iter_dir/grid_search_dmodel/grid_search.out | awk -F"=" '{print \$2}')


# dmodel_dir=$iter_dir/kernel
# dmodel_suffix="_dmodel_scale"

# output_dmodel=1
# out_dir=$iter_dir/model_updated
# mkdir \$out_dir

# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_add_dmodel_alpha_beta_phi_xi_eta \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   $model_dir \$dmodel_dir \$dmodel_suffix \
#   \$step_length \
#   $model_update_min_alpha $model_update_max_alpha \
#   $model_update_min_beta  $model_update_max_beta  \
#   $model_update_min_phi   $model_update_max_phi   \
#   $model_update_min_xi    $model_update_max_xi    \
#   $model_update_min_eta   $model_update_max_eta   \
#   \$output_dmodel \
#   \$out_dir

# # link reference model
# ln -sf $model_REF_dir/proc*_reg1_v[ps]0.bin \$out_dir
# ln -sf $model_REF_dir/proc*_reg1_rho0.bin \$out_dir

# # convert to gll model
# $slurm_mpiexec $sem_utils_dir/bin/xsem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh \
#   $sem_nproc $mesh_dir/DATABASES_MPI \$out_dir  \$out_dir

# # scale rho to beta
# $slurm_mpiexec $sem_utils_dir/bin/xsem_gll_scale_rho_to_beta \
#   $sem_nproc $mesh_dir/DATABASES_MPI \
#   \$out_dir \$out_dir $model_scale_rho_to_beta

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== model_random
# cat <<EOF > $model_random_job
# #!/bin/bash
# #SBATCH -J model_random
# #SBATCH -o ${model_random_job}.o%j
# #SBATCH -N 1
# #SBATCH -n $nproc_per_node
# #SBATCH -p $slurm_partition
# #SBATCH -t $slurm_timelimit_misfit
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# out_dir=model_random
# #rm -rf \$out_dir
# mkdir \$out_dir

# #====== make randomly perturbed model

# min_value=-0.01
# max_value=0.01

# ${slurm_mpiexec} $sem_utils_dir/bin/xsem_gll_random_perturb_model \
#   $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \
#   alpha,beta \$min_value \$max_value \
#   \$out_dir

# #------ link phi,xi,eta,rho from model
# ln -sf $model_dir/proc*_reg1_phi.bin \$out_dir/
# ln -sf $model_dir/proc*_reg1_xi.bin \$out_dir/
# ln -sf $model_dir/proc*_reg1_eta.bin \$out_dir/
# ln -sf $model_dir/proc*_reg1_rho.bin \$out_dir/

# #====== convert to gll model

# # link reference model
# ln -sf $model_REF_dir/proc*_reg1_v[ps]0.bin \$out_dir

# $slurm_mpiexec $sem_utils_dir/bin/xsem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh \
#   $sem_nproc $mesh_dir/DATABASES_MPI \$out_dir  \$out_dir

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== mesh_hess
# cat <<EOF > $mesh_hess_job
# #!/bin/bash
# #SBATCH -J mesh_hess
# #SBATCH -o ${mesh_hess_job}.o%j
# #SBATCH -N $slurm_nnode
# #SBATCH -n $slurm_nproc
# #SBATCH -p $slurm_partition
# #SBATCH -t $slurm_timelimit_mesh
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# for hess_tag in $hess_model_names
# do
#   echo "------ hess model: \${hess_tag}"

#   mesh_dir=$iter_dir/mesh_\${hess_tag}
#   model_dir=$iter_dir/model_\${hess_tag}

#   if [ ! -d "\$model_dir" ]
#   then
#     echo "[ERROR] \$model_dir does not exist!"
#     exit -1
#   fi

#   #rm -rf \$mesh_dir
#   mkdir -p \$mesh_dir

#   cd \$mesh_dir
#   mkdir DATA DATABASES_MPI OUTPUT_FILES

#   cd \$mesh_dir/DATA
#   ln -sf $sem_build_dir/DATA/* \$mesh_dir/DATA/
#   rm Par_file GLL CMTSOLUTION
#   ln -sf \$model_dir GLL
#   cp -L $sem_config_dir/DATA/Par_file .
#   cp -L $sem_config_dir/DATA/CMTSOLUTION .
#   cp -L Par_file CMTSOLUTION \$mesh_dir/OUTPUT_FILES/

#   sed -i "/^MODEL/s/=.*/= GLL/" \$mesh_dir/DATA/Par_file

#   cd \$mesh_dir
#   ${slurm_mpiexec} $sem_build_dir/bin/xmeshfem3D

# done

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# #====== hess_sum
# cat <<EOF > $hess_sum_job
# #!/bin/bash
# #SBATCH -J hess_sum
# #SBATCH -o ${hess_sum_job}.o%j
# #SBATCH -N $slurm_nnode
# #SBATCH -n $slurm_nproc
# #SBATCH -p $slurm_partition
# ##SBATCH -t $slurm_timelimit_misfit
# #SBATCH -t 12:00:00
# #SBATCH --mail-user=kai.tao@utexas.edu
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end

# echo
# echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# #====== process each event kernel
# for event_id in \$(awk -F"|" 'NF&&\$1!~/#/{print \$9}' $event_list)
# do
#   echo "====== link source/receiver mask \$event_id"

#   event_dir=$iter_dir/\$event_id

#   for hess_tag in $hess_model_names
#   do
#     echo "------ hess dmodel: \${hess_tag}"
#     ln -sf \$event_dir/source_receiver_mask/*_mask.bin \$event_dir/output_kernel_\${hess_tag}/kernel
#   done

# done

# echo "====== sum up event kernels [\$(date)]"

# kernel_dir=$iter_dir/kernel

# for hess_tag in $hess_model_names
# do

#   echo ====== \$hess_tag

#   hess_dir=${iter_dir}/hess_\${hess_tag}
#   #rm -rf \$hess_dir
#   mkdir \$hess_dir

#   awk -F"|" 'NF&&\$1!~/#/{printf "%s/%s/output_kernel_%s/kernel\\n", a,\$9,b}' \
#     a="$iter_dir" b="\$hess_tag" $event_list > \$hess_dir/kernel_dir.list

#   echo ------ sum up cijkl,rho_kernel with source and receiver mask

#   ${slurm_mpiexec} $sem_utils_dir/bin/xsem_sum_event_kernels_cijkl \
#     $sem_nproc $mesh_dir/DATABASES_MPI \
#     \$hess_dir/kernel_dir.list cijkl_kernel \
#     1 "mask" \
#     \$hess_dir cijkl_kernel

#   ${slurm_mpiexec} $sem_utils_dir/bin/xsem_sum_event_kernels_1 \
#     $sem_nproc $mesh_dir/DATABASES_MPI \
#     \$hess_dir/kernel_dir.list rho_kernel \
#     1 "mask" \
#     \$hess_dir rho_kernel

#   echo "------ convert cijkl,rho_kernel to aijkl,rhoprime kernel [\$(date)]"

#   ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \
#     $sem_nproc $mesh_dir/DATABASES_MPI ${iter_dir}/model_\${hess_tag} \
#     \$hess_dir \
#     \$hess_dir

#   echo "------ reduce aijkl_kernel to alpha,beta,phi,xi,eta_kernel [\$(date)]"

#   ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_aijkl_to_tiso_in_alpha_beta_phi_xi_eta \
#     $sem_nproc $mesh_dir/DATABASES_MPI ${iter_dir}/model_\${hess_tag} \
#     \$hess_dir \
#     \$hess_dir

#   echo "------ get Hess*dmodel"
#   for kernel_tag in alpha beta phi xi eta rhoprime
#   do
#     #--- Hess*dmodel ~  kernel_dmodel - kernel
#     $slurm_mpiexec $sem_utils_dir/bin/xsem_math \
#       $sem_nproc $mesh_dir/DATABASES_MPI \
#       \$hess_dir \${kernel_tag}_kernel \
#       \$kernel_dir \${kernel_tag}_kernel \
#       "sub" \
#       \$hess_dir \${kernel_tag}_Hdm

#     #--- (Hess*random)^2
#     $slurm_mpiexec $sem_utils_dir/bin/xsem_math \
#       $sem_nproc $mesh_dir/DATABASES_MPI \
#       \$hess_dir \${kernel_tag}_Hdm \
#       \$hess_dir \${kernel_tag}_Hdm \
#       "mul" \
#       \$hess_dir \${kernel_tag}_Hdm_sq
#   done

#   echo ------ precondition kernel
#   for kernel_tag in alpha beta phi xi eta rhoprime
#   do
#     ${slurm_mpiexec} $sem_utils_dir/bin/xsem_math \
#       $sem_nproc $mesh_dir/DATABASES_MPI \
#       \$hess_dir \${kernel_tag}_Hdm \
#       $precond_dir inv_hess_diag \
#       "mul" \
#       \$hess_dir \${kernel_tag}_Hdm_precond
#   done

#   model_tags=alpha_hess_mask_sq,beta_hess_mask_sq,xi_hess_mask_sq,rhoprime_hess_mask_sq

#   $slurm_mpiexec $sem_utils_dir/bin/xsem_smooth \
#     $sem_nproc $mesh_dir/DATABASES_MPI \
#     \$hess_dir \${model_tags} \
#     $hess_smooth_1sigma_h $hess_smooth_1sigma_v \
#     \$hess_dir "_smooth"

#   for kernel_tag in alpha beta xi rhoprime
#   do
#     echo ====== \$kernel_tag

#     $slurm_mpiexec $sem_utils_dir/bin/xsem_math_unary \
#       $sem_nproc $mesh_dir/DATABASES_MPI \
#       \$hess_dir \${kernel_tag}_hess_mask_sq_smooth \
#       "sqrt" \
#       \$hess_dir \${kernel_tag}_hess_mask_sq_smooth_sqrt

#     $slurm_mpiexec $sem_utils_dir/bin/xsem_inverse_hess_diag_water_level \
#       $sem_nproc $mesh_dir/DATABASES_MPI \
#       \$hess_dir \${kernel_tag}_hess_mask_sq_smooth_sqrt \
#       $hess_inverse_nbin $hess_inverse_threshold_percentage \
#       \$hess_dir \${kernel_tag}_inv_hess_mask
#   done

# done

# echo
# echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
# echo

# EOF

# ##END
