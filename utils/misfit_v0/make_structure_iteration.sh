#!/bin/bash
# Make jobs files for slurm 
# structure inversion

# after this put correct CMTSOLUTION into event_id/DATA/ 

wkdir=$(pwd)
sem_utils=/home1/03244/ktao/seiscode/sem_utils
nproc=256

event_id=${1:?[arg]need event_id}

# link the required directories to your wkdir
specfem_dir=$wkdir/specfem3d_globe # bin/x...
mesh_dir=$wkdir/mesh # DATABASES_MPI/proc*_reg1_solver_data.bin
mesh_perturb_dir=$wkdir/mesh_perturb # DATABASES_MPI/proc*_reg1_solver_data.bin
#model_dir=$wkdir/mesh/DATABASES_MPI # proc*_reg1_vph,vpv,vsv,vsh,eta,rho.bin
#!!! model files should reside in mesh_dir/DATABASES_MPI
data_dir=$wkdir/events # <event_id>/data,dis
utils_dir=$wkdir/utils # sem_utils/utils/misfit_v0

# get the full path
specfem_dir=$(readlink -f $specfem_dir)
mesh_dir=$(readlink -f $mesh_dir)
mesh_perturb_dir=$(readlink -f $mesh_perturb_dir)
#model_dir=$(readlink -f $model_dir)
data_dir=$(readlink -f $data_dir)
utils_dir=$(readlink -f $utils_dir)

#====== define variables
# directories
event_dir=$wkdir/$event_id
misfit_dir=$event_dir/misfit
figure_dir=$misfit_dir/figure
slurm_dir=$event_dir/slurm
# job scripts for slurm
mkdir -p $slurm_dir
syn_job=$slurm_dir/syn.job
misfit_job=$slurm_dir/misfit.job
kernel_job=$slurm_dir/kernel.job
hess_job=$slurm_dir/hess.job
precond_job=$slurm_dir/precond.job
perturb_job=$slurm_dir/perturb.job
search_job=$slurm_dir/search.job
# database file
mkdir -p $misfit_dir
db_file=$misfit_dir/misfit.pkl

# misfit par file
misfit_par=$event_dir/DATA/misfit_par.py

#====== syn: forward simulation
cat <<EOF > $syn_job
#!/bin/bash
#SBATCH -J ${event_id}.syn
#SBATCH -o $syn_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 01:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=output_syn

mkdir -p $event_dir/DATA
cd $event_dir/DATA

cp $data_dir/$event_id/data/STATIONS .

cp $mesh_dir/DATA/Par_file .
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .true./" Par_file

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

ibrun $specfem_dir/bin/xspecfem3D

mkdir $event_dir/\$out_dir/sac
mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== misfit
cat <<EOF > $misfit_job
#!/bin/bash
#SBATCH -J ${event_id}.misfit
#SBATCH -o $misfit_job.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -p normal
#SBATCH -t 01:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir

rm -rf $misfit_dir
mkdir -p $misfit_dir
$utils_dir/read_data.py \
  $misfit_par \
  $db_file \
  $event_dir/DATA/CMTSOLUTION \
  $data_dir/$event_id/data/channel.txt \
  $event_dir/output_syn/sac \
  $data_dir/$event_id/dis

$utils_dir/measure_misfit.py $misfit_par $db_file

$utils_dir/output_misfit.py $db_file $misfit_dir/misfit.txt

rm -rf $figure_dir
mkdir -p $figure_dir
$utils_dir/plot_misfit.py $misfit_par $db_file $figure_dir

#------ adjoint source for kernel simulation
rm -rf $event_dir/adj_kernel
mkdir -p $event_dir/adj_kernel
$utils_dir/output_adj.py $db_file $event_dir/adj_kernel

# make STATIONS_ADJOINT
cd $event_dir/adj_kernel
ls *Z.adj | sed 's/..Z\.adj$//' |\
  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
grep -f $event_dir/adj_kernel/grep_pattern $event_dir/DATA/STATIONS \
  > $event_dir/adj_kernel/STATIONS_ADJOINT

#------ adjoint source for hessian simulation
rm -rf $event_dir/adj_hess
mkdir -p $event_dir/adj_hess
$utils_dir/output_adj_hess.py $db_file $event_dir/adj_hess

# make STATIONS_ADJOINT
cd $event_dir/adj_hess
ls *Z.adj | sed 's/..Z\.adj$//' |\
  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
grep -f $event_dir/adj_hess/grep_pattern $event_dir/DATA/STATIONS \
  > $event_dir/adj_hess/STATIONS_ADJOINT

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== kernel simulation
cat <<EOF > $kernel_job
#!/bin/bash
#SBATCH -J ${event_id}.kernel
#SBATCH -o $kernel_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 02:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=output_kernel

cd $event_dir/DATA

sed -i "/^SIMULATION_TYPE/s/=.*/= 3/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^ANISOTROPIC_KL/s/=.*/= .true./" Par_file
sed -i "/^SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
sed -i "/^APPROXIMATE_HESS_KL/s/=.*/= .false./" Par_file

cp $event_dir/adj_kernel/STATIONS_ADJOINT $event_dir/DATA/
rm -rf $event_dir/SEM
ln -s $event_dir/adj_kernel $event_dir/SEM

cd $event_dir

rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

cd $event_dir
ibrun $specfem_dir/bin/xspecfem3D

mkdir $event_dir/\$out_dir/sac
mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

mkdir $event_dir/\$out_dir/kernel
mv $event_dir/DATABASES_MPI/*reg1_cijkl_kernel.bin $event_dir/\$out_dir/kernel/
mv $event_dir/DATABASES_MPI/*reg1_rho_kernel.bin $event_dir/\$out_dir/kernel/

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== hessian simulation
cat <<EOF > $hess_job
#!/bin/bash
#SBATCH -J ${event_id}.hess
#SBATCH -o $hess_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 02:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=output_hess

cd $event_dir/DATA

sed -i "/^SIMULATION_TYPE/s/=.*/= 3/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^ANISOTROPIC_KL/s/=.*/= .true./" Par_file
sed -i "/^SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
sed -i "/^APPROXIMATE_HESS_KL/s/=.*/= .false./" Par_file

cp $event_dir/adj_hess/STATIONS_ADJOINT $event_dir/DATA/
rm -rf $event_dir/SEM
ln -s $event_dir/adj_hess $event_dir/SEM

cd $event_dir

rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

cd $event_dir
ibrun $specfem_dir/bin/xspecfem3D

mkdir $event_dir/\$out_dir/sac
mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

mkdir $event_dir/\$out_dir/kernel
mv $event_dir/DATABASES_MPI/*reg1_cijkl_kernel.bin $event_dir/\$out_dir/kernel/
mv $event_dir/DATABASES_MPI/*reg1_rho_kernel.bin $event_dir/\$out_dir/kernel/

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== kernel preconditioning
cat <<EOF > $precond_job
#!/bin/bash
#SBATCH -J ${event_id}.precond
#SBATCH -o ${precond_job}.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 01:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=$event_dir/kernel_precond
mkdir \$out_dir

echo "====== convert cijkl to aijkl kernel [\$(date)]"
ibrun $sem_utils/bin/xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \
  $nproc $mesh_dir/DATABASES_MPI $mesh_dir/DATABASES_MPI \
  $event_dir/output_kernel/kernel \
  \$out_dir

echo "====== reduce aijkl kernel [\$(date)]"
ibrun $sem_utils/bin/xsem_kernel_aijkl_to_vti_3pars \
  $nproc $mesh_dir/DATABASES_MPI \
  \$out_dir \
  \$out_dir

#echo "====== random kernel to hessian diagonal [\$(date)]"
#ibrun $sem_utils/bin/xsem_hess_diag_sum_random_adjoint_kernel \
#  $nproc $mesh_dir/DATABASES_MPI \
#  $mesh_dir/DATABASES_MPI \
#  $event_dir/output_hess/kernel \
#  \$out_dir

#echo "====== smooth hess diagonal [\$(date)]"
#sigma_h=50
#sigma_v=20
#
#model_tags=sum_hess_diag
#
#ibrun $sem_utils/bin/xsem_smooth \
#  $nproc $mesh_dir/DATABASES_MPI \$out_dir \
#  \$model_tags \$sigma_h \$sigma_v \$out_dir "_smooth"

echo "====== kernel precondition [\$(date)]"
eps=0.001

kernel_tags=vp2_kernel,vsv2_kernel,vsh2_kernel
hess_tag=sum_hess_diag_smooth

ibrun $sem_utils/bin/xsem_kernel_divide_hess_water_level \
  $nproc $mesh_dir/DATABASES_MPI \$out_dir \
  \$kernel_tags \$out_dir \$hess_tag \
  \$eps \$out_dir "_precond"

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== perturb: forward simulation of perturbed model
cat <<EOF > $perturb_job
#!/bin/bash
#SBATCH -J ${event_id}.perturb
#SBATCH -o $perturb_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 02:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir/DATA
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
 
for dmodel in dvp dvsv dvsh
do

  out_dir=output_\${dmodel}
 
  rm -rf $event_dir/DATABASES_MPI
  mkdir $event_dir/DATABASES_MPI
  ln -s $wkdir/mesh_\${dmodel}/DATABASES_MPI/*.bin $event_dir/DATABASES_MPI
  
  cd $event_dir

  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES
  
  cp $wkdir/mesh_\${dmodel}/OUTPUT_FILES/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES
  
  ibrun $specfem_dir/bin/xspecfem3D
  
  mkdir $event_dir/\$out_dir/sac
  mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== search: cc linearized seismograms for chosen step sizes
cat <<EOF > $search_job
#!/bin/bash
#SBATCH -J ${event_id}.search
#SBATCH -o $search_job.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -p normal
#SBATCH -t 01:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

$utils_dir/waveform_der_dmodel.py $misfit_par $db_file $event_dir/output_dvp/sac vp
$utils_dir/waveform_der_dmodel.py $misfit_par $db_file $event_dir/output_dvsv/sac vsv
$utils_dir/waveform_der_dmodel.py $misfit_par $db_file $event_dir/output_dvsh/sac vsh

$utils_dir/grid_search_dvsh.py $misfit_par $db_file $misfit_dir/grid_search_dvsh.txt
$utils_dir/grid_search_dvp_dvsv.py $misfit_par $db_file $misfit_dir/grid_search_dvp_dvsv.txt

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF
