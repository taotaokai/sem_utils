#!/bin/bash
# Make jobs files for slurm 
# structure hessian estimation

wkdir=$(pwd)
sem_utils=/home1/03244/ktao/seiscode/sem_utils
nproc=256

event_id=${1:?[arg]need event_id}
model_name=${2:?[arg]need model_name}

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
#syn_job=$slurm_dir/syn.job
#misfit_job=$slurm_dir/misfit.job
#kernel_job=$slurm_dir/kernel.job
#hess_job=$slurm_dir/hess.job
#precond_job=$slurm_dir/precond.job
perturb_job=$slurm_dir/perturb.job
#search_job=$slurm_dir/search.job
#hess_diag_job=$slurm_dir/hess_diag.job
hess_adj_job=$slurm_dir/hess_adj.job
hess_kernel_job=$slurm_dir/hess_kernel.job
# database file
mkdir -p $misfit_dir
db_file=$misfit_dir/misfit.pkl
# station file
station_file=$event_dir/DATA/STATIONS
if [ ! -f "$station_file" ]
then
  echo "[ERROR] $station_file does NOT exist!"
  exit -1
fi
# cmt file
cmt_file=$event_dir/DATA/CMTSOLUTION.init
if [ ! -f "$cmt_file" ]
then
  echo "[ERROR] $cmt_file does NOT exist!"
  exit -1
fi
# misfit par file
misfit_par=$event_dir/DATA/misfit_par.py
if [ ! -f "$misfit_par" ]
then
  echo "[ERROR] $misfit_par does NOT exist!"
  exit -1
fi

#====== perturb: forward simulation of perturbed model
cat <<EOF > $perturb_job
#!/bin/bash
#SBATCH -J ${event_id}.perturb
#SBATCH -o $perturb_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 00:50:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir/DATA
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
 
for dmodel in $model_name
do

  out_dir=output_perturb_\${dmodel}
 
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

#====== hess_adj
cat <<EOF > $hess_adj_job
#!/bin/bash
#SBATCH -J ${event_id}.hess_adj
#SBATCH -o $hess_adj_job.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -p normal
#SBATCH -t 00:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

$utils_dir/waveform_der_dmodel.py $misfit_par $db_file $event_dir/output_perturb_${model_name}/sac $model_name

rm -rf $event_dir/adj_hess_${model_name}
mkdir -p $event_dir/adj_hess_${model_name}
$utils_dir/output_adj_hess_model_product.py $db_file ${model_name} $event_dir/adj_hess_${model_name}

# make STATIONS_ADJOINT
cd $event_dir/adj_hess_${model_name}
ls *Z.adj | sed 's/..Z\.adj$//' |\
  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
grep -f grep_pattern $event_dir/DATA/STATIONS > STATIONS_ADJOINT

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== hess_kernel simulation
cat <<EOF > $hess_kernel_job
#!/bin/bash
#SBATCH -J ${event_id}.hess_kernel
#SBATCH -o $hess_kernel_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 01:20:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

#model_name=cosine_100km_010

out_dir=output_hess_${model_name}

cd $event_dir/DATA
sed -i "/^SIMULATION_TYPE/s/=.*/= 3/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^ANISOTROPIC_KL/s/=.*/= .true./" Par_file
sed -i "/^SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
sed -i "/^APPROXIMATE_HESS_KL/s/=.*/= .false./" Par_file

cp -f $event_dir/adj_hess_${model_name}/STATIONS_ADJOINT $event_dir/DATA/
rm -rf $event_dir/SEM
ln -s $event_dir/adj_hess_${model_name} $event_dir/SEM

cd $event_dir

rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

rm -rf $event_dir/DATABASES_MPI
mkdir $event_dir/DATABASES_MPI
ln -s $event_dir/forward_saved_frames/*.bin $event_dir/DATABASES_MPI
 
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
