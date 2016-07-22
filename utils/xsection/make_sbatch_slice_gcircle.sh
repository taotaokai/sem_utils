#!/bin/bash

# interpolate SEM gll models to create gcircle slices

#====== command line args
sem_utils=${1:?[arg]need sem_utils dir for bin/xsem_slice_gcircle}
mesh_dir=${2:?[arg]need mesh_dir for proc*_reg1_solver_data.bin}
model_dir=${3:?[arg]need model_dir for proc*_reg1_<model_name>.bin}
slice_list=${4:?[arg]need slice_list}
model_names=${5:?[arg]need model names, e.g. vsv,vsh,vpv,vph,rho,eta}
nc_dir=${6:?[arg]need output directory for .nc files}
job_file=${7:?[arg]need output job_file}
#mpi_exec=${6:?[arg]need mpi_exec, e.g. ibrun or mpirun}

# on TACC:lonestar5
mpi_exec=ibrun

# check input parameters
sem_utils=$(readlink -f $sem_utils)
if [ ! -d "$sem_utils" ]
then
  echo "[ERROR] sem_utils not found: " $sem_utils
  exit -1
fi

mesh_dir=$(readlink -f $mesh_dir)
if [ ! -d "$mesh_dir" ]
then
  echo "[ERROR] mesh_dir not found: " $mesh_dir
  exit -1
fi

model_dir=$(readlink -f $model_dir)
if [ ! -d "$model_dir" ]
then
  echo "[ERROR] model_dir not found: " $model_dir
  exit -1
fi

slice_list=$(readlink -f $slice_list)
if [ ! -f "$slice_list" ]
then
  echo "[ERROR] slice_list not found: " $slice_list
  exit -1
fi

# get some mesh info
nproc=$(ls $mesh_dir/proc*_reg1_solver_data.bin | wc -l)

#====== create job script
cat <<EOF > $job_file
#!/bin/bash
#SBATCH -J slice_gcircle
#SBATCH -o ${job_file}.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 01:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

mkdir $nc_dir

EOF

grep -v "^#" $slice_list |\
while read lat0 lon0 azimuth theta0 theta1 ntheta r0 r1 nr fname
do

cat<<EOF >> $job_file
echo
echo \$(date)
echo "# $fname: $lat0 $lon0 $azimuth $theta0 $theta1 $ntheta $r0 $r1 $nr "
echo
${mpi_exec} \
    $sem_utils/bin/xsem_slice_gcircle \
    $nproc \
    $mesh_dir \
    $model_dir \
    $model_names \
    $lat0 $lon0 $azimuth \
    $theta0 $theta1 $ntheta \
    $r0 $r1 $nr \
    $nc_dir/${fname}.nc

EOF

cat<<EOF >> $job_file
echo
echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

done
