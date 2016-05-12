#!/bin/bash
# Make jobs files for slurm 

wkdir=$(pwd)

event_id=${1:?[arg]need event_id}
iter_num=${2:?[arg]need iter_num}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe

#====== define variables
iter_num=$(echo $iter_num | awk '{printf "%02d",$1}')
iter_prev=$(echo $iter_num | awk '{printf "%02d",$1-1}')
# directories
event_dir=$wkdir/$event_id
misfit_dir=$event_dir/misfit
figure_dir=$misfit_dir/figure.iter${iter_num}
slurm_dir=$event_dir/slurm
# job scripts for slurm
mkdir -p $slurm_dir
green_job=$slurm_dir/iter${iter_num}.green.job
misfit_job=$slurm_dir/iter${iter_num}.misfit.job
srcfrechet_job=$slurm_dir/iter${iter_num}.srcfrechet.job
dgreen_job=$slurm_dir/iter${iter_num}.dgreen.job
search_job=$slurm_dir/iter${iter_num}.search.job
# database file
mkdir -p $misfit_dir
db_file=$misfit_dir/misfit.pkl
# initial cmt file
if [ $iter_num -eq 0 ]
then
  cmt_file=$event_dir/DATA/CMTSOLUTION.harvard.ECEF
else
  cmt_file=$event_dir/DATA/CMTSOLUTION.iter${iter_prev}
fi
if [ ! -f "$cmt_file" ]
then
  echo "[ERROR] $cmt_file does NOT exist!"
  exit -1
fi

#====== green's function
cat <<EOF > $green_job
#!/bin/bash
#SBATCH -J ${event_id}.green.iter$iter_num
#SBATCH -o $green_job.o%j
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

out_dir=output_green

mkdir -p $event_dir/DATA
cd $event_dir/DATA

rm CMTSOLUTION
cp $cmt_file CMTSOLUTION
sed -i "/^tau(s)/s/.*/tau(s):            +0.00000000E+00/" CMTSOLUTION

cp $data_dir/$event_id/data/STATIONS .

cp $mesh_dir/DATA/Par_file .
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file

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

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== misfit
cat <<EOF > $misfit_job
#!/bin/bash
#SBATCH -J ${event_id}.misfit.iter$iter_num
#SBATCH -o $misfit_job.o%j
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

cd $event_dir

mkdir -p $misfit_dir
$utils_dir/read_data.py \
  $db_file \
  $cmt_file \
  $data_dir/$event_id/data/channel.txt \
  $event_dir/output_green \
  $data_dir/$event_id/dis

$utils_dir/measure_misfit.py $db_file

$utils_dir/output_misfit.py $db_file $misfit_dir/misfit.iter${iter_num}.txt

mkdir -p $figure_dir
$utils_dir/plot_misfit.py $db_file $figure_dir

rm -rf $event_dir/SEM
mkdir -p $event_dir/SEM
$utils_dir/output_adj.py $db_file $event_dir/SEM

# make STATIONS_ADJOINT
cd $event_dir/SEM
ls *Z.adj | sed 's/..Z\.adj$//' |\
  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
grep -f $event_dir/SEM/grep_pattern $event_dir/DATA/STATIONS \
  > $event_dir/DATA/STATIONS_ADJOINT

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== derivative of misfit function
cat <<EOF > $srcfrechet_job
#!/bin/bash
#SBATCH -J ${event_id}.srcfrechet.iter$iter_num
#SBATCH -o $srcfrechet_job.o%j
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

out_dir=output_srcfrechet

cd $event_dir/DATA

rm CMTSOLUTION
cp $cmt_file CMTSOLUTION
sed -i "/^tau(s)/s/.*/tau(s):            +0.00000000E+00/" CMTSOLUTION

sed -i "/^SIMULATION_TYPE/s/=.*/= 2/" Par_file

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

mv DATABASES_MPI/*.sem OUTPUT_FILES

echo "make_cmt_der.py [\$(date)]"
$utils_dir/make_cmt_der.py \
  $db_file \
  $event_dir/output_srcfrechet/src_frechet.000001 \
  $event_dir/DATA

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== derivatives of green's function: dxs, dmt
cat <<EOF > $dgreen_job
#!/bin/bash
#SBATCH -J ${event_id}.dgreen.iter$iter_num
#SBATCH -o $dgreen_job.o%j
#SBATCH -N 11
#SBATCH -n 256
#SBATCH -p normal
#SBATCH -t 01:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

for tag in dxs #dmt
do
  echo "====== \$tag"
  out_dir=output_\$tag
  cmt_file=CMTSOLUTION.\$tag

  cd $event_dir/DATA
  rm CMTSOLUTION
  cp -L \$cmt_file CMTSOLUTION
  
  sed -i "/^tau(s)/s/.*/tau(s):            +0.00000000E+00/" CMTSOLUTION 
  sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
  
  cd $event_dir
  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES
  
  cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES
  
  cd $event_dir
  ibrun $specfem_dir/bin/xspecfem3D
done
  
echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== search source parameters
cat <<EOF > $search_job
#!/bin/bash
#SBATCH -J ${event_id}.search.iter$iter_num
#SBATCH -o $search_job.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -p normal
#SBATCH -t 00:40:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir 

# read derivatives of green's fuction 
$utils_dir/waveform_der_xs.py $db_file

# grid search of source model
$utils_dir/search1d_t0_tau_xs.py \
  $db_file DATA/CMTSOLUTION.iter${iter_num} DATA/search.iter${iter_num}.log

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF