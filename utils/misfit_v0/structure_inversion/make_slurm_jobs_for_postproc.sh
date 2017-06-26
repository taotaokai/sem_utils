#!/bin/bash
# Make jobs files for slurm 
# structure inversion

control_file=${1:?[arg]need control_file}
event_list=${2:?[arg]need event_list}

#====== source control_file
source $control_file
event_list=$(readlink -f $event_list)
 
#====== define directories
mesh_dir=$iter_dir/mesh # DATABASES_MPI/proc*_reg1_solver_data.bin
#hess_dir=$iter_dir/hess_sum
mesh_perturb_dir=$iter_dir/mesh_perturb # DATABASES_MPI/proc*_reg1_solver_data.bin
#model_dir=$iter_dir/model # proc*_reg1_vph,vpv,vsv,vsh,eta,rho,dlnvs,kappa,eta,gamma.bin
model_random_dir=$iter_dir/model_random # for model with small random perturbations
kernel_dir=$iter_dir/kernel # DATABASES_MPI/proc*_reg1_solver_data.bin
slurm_dir=$iter_dir/slurm
mkdir -p $slurm_dir
#utils_dir=$sem_utils_dir/utils/misfit_v0 # sem_utils/utils/misfit_v0

#====== define variables
mesh_job=$slurm_dir/mesh.job
kernel_sum_job=$slurm_dir/kernel_sum.job
pcg_dmodel_job=$slurm_dir/pcg_dmodel.job
model_perturb_job=$slurm_dir/model_perturb.job
mesh_perturb_job=$slurm_dir/mesh_perturb.job
model_update_job=$slurm_dir/model_update.job

model_random_job=$slurm_dir/model_random.job
mesh_random_job=$slurm_dir/mesh_random.job
hess_sum_job=$slurm_dir/hess_sum.job

#====== mesh
cat <<EOF > $mesh_job
#!/bin/bash
#SBATCH -J mesh
#SBATCH -o ${mesh_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

mesh_dir="$mesh_dir"
model_dir="$model_dir"

if [ ! -d "\$model_dir" ]
then
  echo "[ERROR] \$model_dir does not exist!"
  exit -1
fi

rm -rf \$mesh_dir
mkdir -p \$mesh_dir

cd \$mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

cd \$mesh_dir/DATA
ln -sf $sem_build_dir/DATA/* \$mesh_dir/DATA/
rm Par_file GLL CMTSOLUTION
ln -sf \$model_dir GLL
cp -L $sem_config_dir/DATA/Par_file .
cp -L $sem_config_dir/DATA/CMTSOLUTION .
cp -L Par_file CMTSOLUTION \$mesh_dir/OUTPUT_FILES/

sed -i "/^MODEL/s/=.*/= GLL/" \$mesh_dir/DATA/Par_file

cd \$mesh_dir
${slurm_mpiexec} $sem_build_dir/bin/xmeshfem3D

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== kernel_sum
cat <<EOF > $kernel_sum_job
#!/bin/bash
#SBATCH -J kernel_sum
#SBATCH -o ${kernel_sum_job}.o%j
#SBATCH -N 1
#SBATCH -n $nproc_per_node
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

#====== process each event kernel
for event_id in \$(awk -F"|" 'NF&&\$1!~/#/{print \$9}' $event_list)
do
  echo "====== process \$event_id"
  event_dir=$iter_dir/\$event_id
  out_dir=\$event_dir/kernel
  mkdir \$out_dir

  echo "------ convert cijkl to aijkl kernel [\$(date)]"
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \
    $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \
    \$event_dir/output_kernel/kernel \
    \$out_dir

  echo "------ reduce aijkl kernel to alpha,beta,xi kernel [\$(date)]"
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_aijkl_to_tiso_in_alpha_beta_xi_scale_phi_eta_to_xi \
    $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \$out_dir $model_scale_phi_to_xi $model_scale_eta_to_xi \$out_dir

  echo "------ make kernel mask [\$(date)]"
  awk 'NR==6{print \$0, a}' a=$source_mask_1sigma_km \$event_dir/output_kernel/source.vtk > \$out_dir/source.xyz
  awk 'NR>=6&&NF==3{print \$0, a}' a=$receiver_mask_1sigma_km \$event_dir/output_kernel/receiver.vtk >> \$out_dir/source.xyz
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_make_gaussian_mask \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    \$out_dir/source.xyz \
    \$out_dir "mask"
done

#====== sum up event kernels
echo "====== sum up event kernels [\$(date)]"
awk -F"|" 'NF&&\$1!~/#/{printf "%s/%s/kernel\\n", a,\$9}' \
  a="$iter_dir" $event_list > kernel_dir.list

out_dir=$iter_dir/kernel_sum
mkdir \$out_dir

# raw kernel
for kernel_tag in alpha beta xi rhoprime
do
  echo ====== raw kernel: \$kernel_tag
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_sum_event_kernels_1 \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    kernel_dir.list \${kernel_tag}_kernel \
    0 "xxx" \
    \$out_dir \${kernel_tag}_kernel
done

# preconditioned kernel
for kernel_tag in alpha beta xi
do
  echo ====== preconditioned kernel: \$kernel_tag

  echo ------ mask source/receiver
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_sum_event_kernels_1 \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    kernel_dir.list \${kernel_tag}_kernel \
    1 "mask" \
    \$out_dir \${kernel_tag}_kernel_mask

  echo ------ inverse hessian diagonals
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_math \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    \$out_dir \${kernel_tag}_kernel_mask \
    $hess_dir \${kernel_tag}_inv_hess_diag \
    "mul" \
    \$out_dir \${kernel_tag}_kernel_precond
done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== pcg_dmodel
cat <<EOF > $pcg_dmodel_job
#!/bin/bash
#SBATCH -J pcg_dmodel
#SBATCH -o ${pcg_dmodel_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

#------ preconditioned conjugate gradient
if [ x$iter_num != x00 ]
then

  current_precond_kernel_names=dlnvs_kernel_precond,kappa_kernel_precond,eps_kernel_precond,gamma_kernel_precond
  current_kernel_names=dlnvs_kernel,kappa_kernel,eps_kernel,gamma_kernel
  previous_kernel_names=dlnvs_kernel,kappa_kernel,eps_kernel,gamma_kernel
  previous_dmodel_names=dlnvs_dmodel,kappa_dmodel,eps_dmodel,gamma_dmodel
  out_dmodel_names=dlnvs_dmodel,kappa_dmodel,eps_dmodel,gamma_dmodel
  
  previous_kernel_dir=$prev_iter_dir/kernel
  previous_dmodel_dir=$prev_iter_dir/model_updated
  
  $slurm_mpiexec $sem_utils_dir/bin/xsem_pcg_dmodel_n \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $kernel_dir  \$current_precond_kernel_names \
    $kernel_dir  \$current_kernel_names \
    \$previous_kernel_dir \$previous_kernel_names \
    \$previous_dmodel_dir \$previous_dmodel_names \
    $cg_type \
    $kernel_dir \$out_dmodel_names

else

  cd $kernel_dir
  ls *_precond.bin | awk -F"_" '{printf "ln -s %s %s_%s_%s_dmodel.bin\\n",\$0,\$1,\$2,\$3}' > ln_sh
  bash ln_sh

fi

#------ scale dmodel
kernel_suffix="_dmodel"
out_suffix="_dmodel_scale"

$slurm_mpiexec $sem_utils_dir/bin/xsem_scale_dmodel_dlnvs_kappa_thomsen_elliptic \
  $sem_nproc $mesh_dir/DATABASES_MPI \
  $kernel_dir \$kernel_suffix \
  $dmodel_max_dlnvs \
  $kernel_dir \$out_suffix

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== model_perturb
cat <<EOF > $model_perturb_job
#!/bin/bash
#SBATCH -J model_perturb
#SBATCH -o ${model_perturb_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=$iter_dir/model_perturb
mkdir \$out_dir

#====== add dmodel
model_names=dlnvs,kappa,eps,gamma
dmodel_names=dlnvs_dmodel_scale,kappa_dmodel_scale,eps_dmodel_scale,gamma_dmodel_scale
${slurm_mpiexec} $sem_utils_dir/bin/xsem_math \
  $sem_nproc $mesh_dir/DATABASES_MPI \
  $model_dir \$model_names \
  $kernel_dir \$dmodel_names \
  "add" \
  \$out_dir \$model_names

#====== convert to gll model
ls $mesh_REF_dir/DATABASES_MPI/proc*_reg1_vsv.bin | awk -F"/" '{a=\$NF;gsub(/\.bin/,"_ref.bin",a); printf "ln -s %s %s/%s\\n",\$0,outdir,a}' outdir=\$out_dir > \$out_dir/link_vsv_ref_sh
bash \$out_dir/link_vsv_ref_sh

${slurm_mpiexec} $sem_utils_dir/bin/xsem_model_dlnvs_kappa_thomsen_elliptic_to_gll \
  $sem_nproc $mesh_dir/DATABASES_MPI \$out_dir x \$out_dir

#====== use reference rho model
ln -s $mesh_REF_dir/DATABASES_MPI/proc*_reg1_rho.bin \$out_dir/

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== mesh_perturb
cat <<EOF > $mesh_perturb_job
#!/bin/bash
#SBATCH -J mesh_perturb
#SBATCH -o ${mesh_perturb_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

mesh_dir=$iter_dir/mesh_perturb
model_dir=$iter_dir/model_perturb

if [ ! -d "\$model_dir" ]
then
  echo "[ERROR] \$model_dir does not exist!"
  exit -1
fi

rm -rf \$mesh_dir
mkdir -p \$mesh_dir

cd \$mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

cd \$mesh_dir/DATA
ln -sf $sem_build_dir/DATA/* \$mesh_dir/DATA/
rm Par_file GLL CMTSOLUTION
ln -sf \$model_dir GLL
cp -L $sem_config_dir/DATA/Par_file .
cp -L $sem_config_dir/DATA/CMTSOLUTION .
cp -L Par_file CMTSOLUTION \$mesh_dir/OUTPUT_FILES/

sed -i "/^MODEL/s/=.*/= GLL/" \$mesh_dir/DATA/Par_file

${slurm_mpiexec} $sem_build_dir/bin/xmeshfem3D

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== model_update
cat <<EOF > $model_update_job
#!/bin/bash
#SBATCH -J model_update
#SBATCH -o ${model_update_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $iter_dir
mkdir grid_search_dmodel
ls */misfit/grid_search_dmodel.txt | awk -F"/" '{printf "cp %s grid_search_dmodel/%s.txt\\n",\$0,\$1}' | bash

cd $iter_dir/grid_search_dmodel
ls *.txt > list
$sem_utils_dir/utils/misfit_v0/plot_grid_search_1d.py list > grid_search.out
step_length=\$(awk -F"=" '{print \$2}' grid_search.out)


dmodel_dir=$iter_dir/kernel_sum
dmodel_suffix="_dmodel_scale"

output_dmodel=1
out_dir=$iter_dir/model_updated
mkdir \$out_dir

${slurm_mpiexec} $sem_utils_dir/bin/xsem_add_dmodel_dlnvs_kappa_thomsen_elliptic \
  $sem_nproc $mesh_dir/DATABASES_MPI \
  $model_dir \$dmodel_dir \$dmodel_suffix \
  \$step_length \
  \$model_update_min_dlnvs \$model_update_max_dlnvs \
  \$model_update_min_kappa \$model_update_max_kappa \
  \$model_update_min_eps   \$model_update_max_eps \
  \$model_update_min_gamma \$model_update_max_gamma \
  \$output_dmodel \
  \$out_dir

# convert to gll model
ls $mesh_REF_dir/DATABASES_MPI/proc*_reg1_vsv.bin | awk -F"/" '{a=\$NF;gsub(/\.bin/,"_ref.bin",a); printf "ln -s %s %s/%s\\n",\$0,outdir,a}' outdir=\$out_dir > \$out_dir/link_vsv_ref_sh
bash \$out_dir/link_vsv_ref_sh

${slurm_mpiexec} $sem_utils_dir/bin/xsem_model_dlnvs_kappa_thomsen_elliptic_to_gll \
  $sem_nproc $mesh_dir/DATABASES_MPI \$out_dir x \$out_dir

# use reference rho model
ln -s $mesh_REF_dir/DATABASES_MPI/proc*_reg1_rho.bin \$out_dir/

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== model_random
cat <<EOF > $model_random_job
#!/bin/bash
#SBATCH -J model_random
#SBATCH -o ${model_random_job}.o%j
#SBATCH -N 1
#SBATCH -n $nproc_per_node
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=model_random
rm -rf \$out_dir
mkdir \$out_dir

#====== make randomly perturbed model

min_value=-0.01
max_value=0.01

${slurm_mpiexec} $sem_utils_dir/bin/xsem_gll_random_perturb_model \
  $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \
  alpha,beta,xi \$min_value \$max_value \
  \$out_dir

#====== convert to gll model

# link reference model
cat /dev/null > \$out_dir/link_sh

ls $mesh_REF_dir/DATABASES_MPI/proc*_reg1_vpv.bin | awk -F"/" '{a=\$NF;gsub(/v\.bin/,"0.bin",a); printf "ln -s %s %s/%s\\n",\$0,outdir,a}' outdir=\$out_dir >> \$out_dir/link_sh

ls $mesh_REF_dir/DATABASES_MPI/proc*_reg1_vsv.bin | awk -F"/" '{a=\$NF;gsub(/v\.bin/,"0.bin",a); printf "ln -s %s %s/%s\\n",\$0,outdir,a}' outdir=\$out_dir >> \$out_dir/link_sh

ls $mesh_REF_dir/DATABASES_MPI/proc*_reg1_rho.bin | awk -F"/" '{a=\$NF;gsub(/\.bin/,"0.bin",a); printf "ln -s %s %s/%s\\n",\$0,outdir,a}' outdir=\$out_dir >> \$out_dir/link_sh

bash \$out_dir/link_sh

$slurm_mpiexec $sem_utils_dir/bin/xsem_gll_alpha_beta_xi_to_tiso_scale_phi_eta_rho_to_xi \
  $sem_nproc $mesh_dir/DATABASES_MPI \$out_dir $model_scale_phi_to_xi $model_scale_eta_to_xi $model_scale_rho_to_beta \$out_dir

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== mesh_random
cat <<EOF > $mesh_random_job
#!/bin/bash
#SBATCH -J mesh_random
#SBATCH -o ${mesh_random_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
#SBATCH -t $slurm_timelimit_misfit
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

mesh_dir=$iter_dir/mesh_random
model_dir=$iter_dir/model_random

if [ ! -d "\$model_dir" ]
then
  echo "[ERROR] \$model_dir does not exist!"
  exit -1
fi

rm -rf \$mesh_dir
mkdir -p \$mesh_dir

cd \$mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

cd \$mesh_dir/DATA
ln -sf $sem_build_dir/DATA/* \$mesh_dir/DATA/
rm Par_file GLL CMTSOLUTION
ln -sf \$model_dir GLL
cp -L $sem_config_dir/DATA/Par_file .
cp -L $sem_config_dir/DATA/CMTSOLUTION .
cp -L Par_file CMTSOLUTION \$mesh_dir/OUTPUT_FILES/

sed -i "/^MODEL/s/=.*/= GLL/" \$mesh_dir/DATA/Par_file

cd \$mesh_dir
${slurm_mpiexec} $sem_build_dir/bin/xmeshfem3D

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== hess_sum
cat <<EOF > $hess_sum_job
#!/bin/bash
#SBATCH -J hess_sum
#SBATCH -o ${hess_sum_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition
##SBATCH -t $slurm_timelimit_misfit
#SBATCH -t 12:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

#====== process each event kernel
for event_id in \$(awk -F"|" 'NF&&\$1!~/#/{print \$9}' $event_list)
do
  echo "====== process \$event_id"

  event_dir=$iter_dir/\$event_id

  #--- kernel
  out_dir=\$event_dir/kernel
  rm -rf \$out_dir
  mkdir \$out_dir

  echo "------ convert cijkl to aijkl kernel [\$(date)]"
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \
    $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \
    \$event_dir/output_kernel/kernel \
    \$out_dir

  echo "------ reduce aijkl kernel to alpha,beta,xi kernel [\$(date)]"
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_aijkl_to_tiso_in_alpha_beta_xi_scale_phi_eta_to_xi \
    $sem_nproc $mesh_dir/DATABASES_MPI $model_dir \$out_dir $model_scale_phi_to_xi $model_scale_eta_to_xi \$out_dir

  echo "------ make kernel mask [\$(date)]"
  awk 'NR==6{print \$0, a}' a=$source_mask_1sigma_km \$event_dir/output_kernel/source.vtk > \$out_dir/source.xyz
  awk 'NR>=6&&NF==3{print \$0, a}' a=$receiver_mask_1sigma_km \$event_dir/output_kernel/receiver.vtk >> \$out_dir/source.xyz
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_make_gaussian_mask \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    \$out_dir/source.xyz \
    \$out_dir "mask"

  #--- kernel_random
  out_dir=\$event_dir/kernel_random
  rm -rf \$out_dir
  mkdir \$out_dir

  echo "------ convert cijkl to aijkl kernel [\$(date)]"
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \
    $sem_nproc $mesh_dir/DATABASES_MPI $model_random_dir \
    \$event_dir/output_kernel_random/kernel \
    \$out_dir

  echo "------ reduce aijkl kernel to alpha,beta,xi kernel [\$(date)]"
  ${slurm_mpiexec} $sem_utils_dir/bin/xsem_kernel_aijkl_to_tiso_in_alpha_beta_xi_scale_phi_eta_to_xi \
    $sem_nproc $mesh_dir/DATABASES_MPI $model_random_dir \$out_dir $model_scale_phi_to_xi $model_scale_eta_to_xi \$out_dir

  ln -sf \$event_dir/kernel/*_mask.bin \$event_dir/kernel_random/

done

#====== sum up event kernels
echo "====== sum up event kernels [\$(date)]"
awk -F"|" 'NF&&\$1!~/#/{printf "%s/%s/kernel\\n", a,\$9}' \
  a="$iter_dir" $event_list > $iter_dir/kernel_dir.list

awk -F"|" 'NF&&\$1!~/#/{printf "%s/%s/kernel_random\\n", a,\$9}' \
  a="$iter_dir" $event_list > $iter_dir/kernel_random_dir.list

rm -rf $hess_dir
mkdir $hess_dir

for kernel_tag in alpha beta xi rhoprime
do
  echo ====== \$kernel_tag

  $slurm_mpiexec $sem_utils_dir/bin/xsem_sum_event_kernels_1 \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $iter_dir/kernel_dir.list \${kernel_tag}_kernel \
    1 "mask" \
    $hess_dir \${kernel_tag}_kernel_mask

  $slurm_mpiexec $sem_utils_dir/bin/xsem_sum_event_kernels_1 \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $iter_dir/kernel_random_dir.list \${kernel_tag}_kernel \
    1 "mask" \
    $hess_dir \${kernel_tag}_kernel_random_mask

  #--- Hess*random ~  kernel_random - kernel
  $slurm_mpiexec $sem_utils_dir/bin/xsem_math \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $hess_dir \${kernel_tag}_kernel_mask \
    $hess_dir \${kernel_tag}_kernel_random_mask \
    "sub" \
    $hess_dir \${kernel_tag}_hess_random_mask

done

# approximated Hessian diagonals
for kernel_tag in alpha beta xi rhoprime
do
  echo ====== \$kernel_tag

  $slurm_mpiexec $sem_utils_dir/bin/xsem_math \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $hess_dir \${kernel_tag}_hess_random_mask \
    $hess_dir \${kernel_tag}_hess_random_mask \
    "mul" \
    $hess_dir \${kernel_tag}_hess_random_mask_diag_sq

done

#model_tags=alpha_hess_random_diag_sq,beta_hess_random_diag_sq,xi_hess_random_diag_sq,rhoprime_hess_random_diag_sq
model_tags=alpha_hess_random_mask_diag_sq,beta_hess_random_mask_diag_sq,xi_hess_random_mask_diag_sq,rhoprime_hess_random_mask_diag_sq

$slurm_mpiexec $sem_utils_dir/bin/xsem_smooth \
  $sem_nproc $mesh_dir/DATABASES_MPI \
  $hess_dir \${model_tags} \
  $hess_smooth_1sigma_h $hess_smooth_1sigma_v \
  $hess_dir "_smooth"

for kernel_tag in alpha beta xi rhoprime
do
  echo ====== \$kernel_tag

  $slurm_mpiexec $sem_utils_dir/bin/xsem_math_unary \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $hess_dir \${kernel_tag}_hess_random_mask_diag_sq_smooth \
    "sqrt" \
    $hess_dir \${kernel_tag}_hess_mask_diag

  $slurm_mpiexec $sem_utils_dir/bin/xsem_inverse_hess_diag_water_level \
    $sem_nproc $mesh_dir/DATABASES_MPI \
    $hess_dir \${kernel_tag}_hess_mask_diag \
    $hess_inverse_nbin $hess_inverse_threshold_percentage \
    $hess_dir \${kernel_tag}_inv_hess_mask_diag

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

##END
