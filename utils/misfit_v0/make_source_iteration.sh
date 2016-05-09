#!/bin/bash
# Make jobs files for slurm 

wkdir=$(pwd)

event_id=${1:?[arg]need event_id}
iter_num=${2:?[arg]need iter_num}

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe

event_dir=$wkdir/$event_id
iter_num=$(echo $iter_num | awk '{printf "%02d",$1}')
iter_prev=$(echo $iter_num | awk '{printf "%02d",$1-1}')

# job scripts
green_job=$event_dir/green.job
misfit_job=$event_dir/misfit.job
srcfrechet_job=$event_dir/srcfrechet.job
dxs_dmt_job=$event_dir/dxs_dmt.job
search_job=$event_dir/search.job

# database file
db_file=$event_dir/misfit.pkl

# check event dir
if [ ! -d "$event_dir" ]
then
  echo "[ERROR] $event_dir does NOT exist!"
  exit -1
fi

# initial cmt file
if [ $iter_num -eq 0 ]
then
  cmt_file=$event_dir/DATA/CMTSOLUTION.harvard.ECEF
else
  cmt_file=$event_dir/DATA/CMTSOLUTION.${iter_prev}
fi
if [ ! -f "$cmt_file" ]
then
  echo "[ERROR] $cmt_file does NOT exist!"
  exit -1
fi

#------ green function
cat <<EOF > $green_job
#!/bin/bash
#SBATCH -J $event_id
#SBATCH -o $green_job.${iter_num}.o%j
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

cd $event_dir/DATA

rm CMTSOLUTION
cp $cmt_file CMTSOLUTION
sed -i "/^tau(s)/s/.*/tau(s):            +0.00000000E+00/" CMTSOLUTION

cp $data_dir/$event_id/data/STATIONS .

cp $mesh_dir/DATA/Par_file.
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file

cd $event_dir
rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

cp DATABASES_MPI/addressing.txt OUTPUT_FILES
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
#SBATCH -J $event_id
#SBATCH -o $misfit_job.${iter_num}.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 00:50:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir

$utils_dir/read_data.py \
  $db_file \
  $cmt_file \
  $data_dir/$event_id/data/channel.txt \
  $event_dir/output_green \
  $data_dir/$event_id/dis

$utils_dir/measure_misfit.py $db_file

rm -rf $event_dir/SEM
mkdir -p $event_dir/SEM
$utils_dir/output_adj.py $db_file $event_dir/SEM

# STATIONS_ADJOINT


echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== derivative of misfit function
cat <<EOF > $srcfrechet_job
#!/bin/bash
#SBATCH -J $event_id
#SBATCH -o $srcfrechet_job.${iter_num}o%j
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
rm -rf \$out_dir OUTPUT_FILES SEM
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

cp DATABASES_MPI/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

ibrun $specfem_dir/bin/xspecfem3D

mv DATABASES_MPI/*.sem OUTPUT_FILES

$utils_dir/make_cmt_der.py \
  $event_dir/misfit.pkl \
  $event_dir/output_srcfrechet/src_frechet.000001 \
  $event_dir/DATA \
  > $event_dir/make_cmt_der.log 2>&1

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== derivatives of green's function: dxs, dmt
cat <<EOF > $event_dir/dxs_dmt.job
#!/bin/bash
#SBATCH -J $event_id
#SBATCH -o $event_dir/dxs_dmt.job.o%j
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

for tag in dxs dmt
do
  echo "====== \$tag"
  out_dir=output_\$tag
  cmt_file=CMTSOLUTION.\$tag

  cd $event_dir/DATA
  rm CMTSOLUTION
  ln -sf \$cmt_file CMTSOLUTION
  
  sed -i "/^tau(s)/s/.*/tau(s):            +0.00000000E+00/" \$cmt_file 
  sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
  
  cd $event_dir
  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES
  
  cp DATABASES_MPI/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES
  
  ibrun $specfem_dir/bin/xspecfem3D
done
  
echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF

#====== grid search
cat <<EOF > $event_dir/search.job
#!/bin/bash
#SBATCH -J $event_id
#SBATCH -o $event_dir/search.job.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 00:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

# read derivatives of green's fuction 
$utils_dir/waveform_der.py $misfit_file

# grid search of source model
$utils_dir/search1d.py $misfit_file

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF