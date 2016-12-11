#!/bin/bash
# Make jobs files for slurm 
# structure inversion

# after this put correct CMTSOLUTION into event_id/DATA/ 

wkdir=$(pwd)

event_id=${1:?[arg]need event_id}

# link the required directories to your wkdir
specfem_dir=$wkdir/specfem3d_globe
mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils

# get the full path
specfem_dir=$(readlink -f $specfem_dir)
mesh_dir=$(readlink -f $mesh_dir)
data_dir=$(readlink -f $data_dir)
utils_dir=$(readlink -f $utils_dir)

#====== define variables
# directories
event_dir=$wkdir/$event_id
misfit_dir=$event_dir/misfit2
figure_dir=$misfit_dir/figure
slurm_dir=$event_dir/slurm
# job scripts for slurm
mkdir -p $slurm_dir
misfit_job=$slurm_dir/misfit2.job
# database file
mkdir -p $misfit_dir
db_file=$misfit_dir/misfit.pkl

# misfit par file
misfit_par=$event_dir/DATA/misfit_par2.py

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
  $db_file \
  $event_dir/DATA/CMTSOLUTION \
  $data_dir/$event_id/data/channel.txt \
  $event_dir/output_syn/sac \
  $data_dir/$event_id/dis

# backup misfit_par file
cp $misfit_par $misfit_dir

$utils_dir/measure_misfit.py $misfit_par $db_file

$utils_dir/output_misfit.py $db_file $misfit_dir/misfit.txt

rm -rf $figure_dir
mkdir -p $figure_dir
$utils_dir/plot_misfit.py $misfit_par $db_file $figure_dir

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF
