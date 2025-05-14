#!/bin/bash

# Make jobs files for slurm
# source inversion

control_file=${1:?[arg]need control_file}
event_id=${2:?[arg]need event_id}

# source control_file
source $control_file

#====== define variables
# directories
event_dir=$iter_dir/$event_id
#mesh_dir=$iter_dir/mesh
misfit_dir=$event_dir/misfit
figure_dir=$misfit_dir/figure
slurm_dir=$event_dir/slurm
# job scripts for slurm
mkdir -p $slurm_dir
green_job=$slurm_dir/green.job
misfit_job=$slurm_dir/misfit.job
srcfrechet_job=$slurm_dir/srcfrechet.job
dgreen_job=$slurm_dir/dgreen.job
search_job=$slurm_dir/search.job

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
cmt_file=$event_dir/DATA/CMTSOLUTION.init
if [ ! -f "$cmt_file" ]
then
  echo "[WARN] $cmt_file does NOT exist!"
  exit -1
fi
# misfit par file
misfit_par=$event_dir/DATA/misfit.yaml
if [ ! -f "$misfit_par" ]
then
  echo "[ERROR] $misfit_par does NOT exist!"
  exit -1
fi

#====== green's function
cat <<EOF > $green_job
#!/bin/bash
#SBATCH -J ${event_id}.green
#SBATCH -o ${green_job}.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition_dcu
#SBATCH -t $slurm_timelimit_forward
#SBATCH $slurm_dcu_extra_args

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=output_green

chmod u+w -R $event_dir/DATA
mkdir -p $event_dir/DATA
cd $event_dir/DATA

rm CMTSOLUTION
cp $cmt_file CMTSOLUTION
sed -i "/^tau(s)/s/.*/tau(s):            +0.0E+00/" CMTSOLUTION

# cp $data_dir/$event_id/STATIONS .

cp $mesh_dir/DATA/Par_file .
sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^USE_ECEF_COORDINATE/s/=.*/= .true./" Par_file
sed -i "/^USE_FORCE_POINT_SOURCE/s/=.*/= .false./" Par_file

rm -rf $event_dir/DATABASES_MPI
mkdir $event_dir/DATABASES_MPI
ln -s $mesh_dir/DATABASES_MPI/*.bin $event_dir/DATABASES_MPI

cd $event_dir
chmod u+w -R \$out_dir
rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

# cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

$slurm_mpiexec $sem_build_dir/bin/xspecfem3D

mkdir $event_dir/\$out_dir/sac
mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

# modify CMTSOLUTION.init with the source location actually used in the simulation
tmpfile=\$(mktemp)
grep -A4 "position of the source that will be used:" $event_dir/\$out_dir/output_solver.txt > \$tmpfile
if [ \$? -ne 0 ]
then
  echo "[ERROR] check if green.job finished OK!"
  exit -1
else
  cp $cmt_file ${cmt_file}.orig
  x=\$(grep "x(m)" \$tmpfile | awk '{printf "%+15.8E", \$2}')
  y=\$(grep "y(m)" \$tmpfile | awk '{printf "%+15.8E", \$2}')
  z=\$(grep "z(m)" \$tmpfile | awk '{printf "%+15.8E", \$2}')
  sed -i "s/x(m).*/x(m):              \$x/"  $cmt_file
  sed -i "s/y(m).*/y(m):              \$y/"  $cmt_file
  sed -i "s/z(m).*/z(m):              \$z/"  $cmt_file
fi
rm \$tmpfile

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
#SBATCH -p $slurm_partition_cpu
#SBATCH -t $slurm_timelimit_misfit

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir

chmod u+w -R $misfit_dir
rm -rf $misfit_dir
mkdir -p $misfit_dir

chmod u+w $event_dir/SEM
rm -rf $event_dir/SEM
mkdir -p $event_dir/SEM

$python_exec $sem_utils_dir/misfit/measure_adj.py \\
  $db_file \\
  $misfit_par \\
  $event_dir/DATA/Par_file \\
  $cmt_file \\
  $data_dir/$event_id/channel.txt \\
  $data_dir/$event_id/data.h5 \\
  $event_dir/output_green/sac \\
  $event_dir/SEM \\
  --cmt_in_ECEF \\
  --syn_is_grn

# make STATIONS_ADJOINT
cd $event_dir/SEM
ls *Z.adj | sed 's/..Z\.adj$//' |\
  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",\$1,\$2,\$3}' > grep_pattern
grep -f $event_dir/SEM/grep_pattern $event_dir/DATA/STATIONS \
  > $event_dir/SEM/STATIONS_ADJOINT

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== source frechet simulation
cat <<EOF > $srcfrechet_job
#!/bin/bash
#SBATCH -J ${event_id}.srcfrechet
#SBATCH -o $srcfrechet_job.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition_dcu
#SBATCH -t $slurm_timelimit_forward
#SBATCH $slurm_dcu_extra_args

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

out_dir=output_srcfrechet

cd $event_dir/DATA

rm CMTSOLUTION
cp $cmt_file CMTSOLUTION
sed -i "/^SIMULATION_TYPE/s/=.*/= 2/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^USE_ECEF_COORDINATE/s/=.*/= .true./" Par_file
sed -i "/^USE_FORCE_POINT_SOURCE/s/=.*/= .false./" Par_file

cd $event_dir

rm -rf \$out_dir OUTPUT_FILES
mkdir \$out_dir
ln -sf \$out_dir OUTPUT_FILES

chmod u+w $event_dir/DATA/STATIONS_ADJOINT
cp SEM/STATIONS_ADJOINT DATA/

cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES

cd $event_dir
$slurm_mpiexec $sem_build_dir/bin/xspecfem3D

# mv DATABASES_MPI/*.sem OUTPUT_FILES

cp $event_dir/\$out_dir/src_frechet.000001 $misfit_dir/srcfrechet

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== partial derivative of green's function
cat <<EOF > $dgreen_job
#!/bin/bash
#SBATCH -J ${event_id}.dgreen
#SBATCH -o $dgreen_job.o%j
#SBATCH -N $slurm_nnode
#SBATCH -n $slurm_nproc
#SBATCH -p $slurm_partition_dcu
#SBATCH -t $slurm_timelimit_dgreen
#SBATCH $slurm_dcu_extra_args

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

# make perturbed CMTSOLUTION
$python_exec $sem_utils_dir/misfit/make_perturbed_cmtsolution.py \\
  $db_file \\
  $misfit_dir/srcfrechet \\
  $misfit_dir/diff_CMTSOLUTION \\
  $misfit_dir/CMTSOLUTION_dxs \\
  $misfit_dir/CMTSOLUTION_dmt

for tag in dxs dmt
do
  echo "====== \$tag"
  out_dir=output_\$tag
  dcmt_file=$misfit_dir/CMTSOLUTION_\$tag

  cd $event_dir/DATA

  rm CMTSOLUTION
  cp -L \$dcmt_file CMTSOLUTION
  sed -i "/^tau(s)/s/.*/tau(s):            +0.00000000E+00/" CMTSOLUTION

  sed -i "/^SIMULATION_TYPE/s/=.*/= 1/" Par_file
  sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
  sed -i "/^USE_ECEF_COORDINATE/s/=.*/= .true./" Par_file
  sed -i "/^USE_FORCE_POINT_SOURCE/s/=.*/= .false./" Par_file

  cd $event_dir
  rm -rf \$out_dir OUTPUT_FILES
  mkdir \$out_dir
  ln -sf \$out_dir OUTPUT_FILES

  # cp $mesh_dir/OUTPUT_FILES/addressing.txt OUTPUT_FILES
  cp -L DATA/Par_file OUTPUT_FILES
  cp -L DATA/STATIONS OUTPUT_FILES
  cp -L DATA/CMTSOLUTION OUTPUT_FILES

  cd $event_dir
  $slurm_mpiexec $sem_build_dir/bin/xspecfem3D

  mkdir $event_dir/\$out_dir/sac
  mv $event_dir/\$out_dir/*.sac $event_dir/\$out_dir/sac

  # modify diff_CMTSOLUTION with the perturbation in source location actually used
  if [ x"\$tag" == x"dxs" ]
  then
    tmpfile=\$(mktemp)
    grep -A4 "position of the source that will be used:" $event_dir/\$out_dir/output_solver.txt > \$tmpfile
    if [ \$? -ne 0 ]
    then
      echo "[ERROR] check if dgreen.job finished OK!"
      exit -1
    fi

    x1=\$(grep "x(m)" \$tmpfile | awk '{printf "%+15.8E", \$2}')
    y1=\$(grep "y(m)" \$tmpfile | awk '{printf "%+15.8E", \$2}')
    z1=\$(grep "z(m)" \$tmpfile | awk '{printf "%+15.8E", \$2}')

    sed -i "s/x(m).*/x(m):              \$x1/"  \$dcmt_file
    sed -i "s/y(m).*/y(m):              \$y1/"  \$dcmt_file
    sed -i "s/z(m).*/z(m):              \$z1/"  \$dcmt_file

    # modify dcmt file
    x0=\$(grep "x(m)" $cmt_file | awk '{printf "%+15.8E", \$2}')
    y0=\$(grep "y(m)" $cmt_file | awk '{printf "%+15.8E", \$2}')
    z0=\$(grep "z(m)" $cmt_file | awk '{printf "%+15.8E", \$2}')

    dx=\$(echo \$x1 \$x0 | awk '{printf "%+15.8E", \$1-\$2}')
    dy=\$(echo \$y1 \$y0 | awk '{printf "%+15.8E", \$1-\$2}')
    dz=\$(echo \$z1 \$z0 | awk '{printf "%+15.8E", \$1-\$2}')

    # original dxs
    dx0=\$(grep "dx(m)" $misfit_dir/diff_CMTSOLUTION | awk '{printf "%+15.8E", \$2}')
    dy0=\$(grep "dy(m)" $misfit_dir/diff_CMTSOLUTION | awk '{printf "%+15.8E", \$2}')
    dz0=\$(grep "dz(m)" $misfit_dir/diff_CMTSOLUTION | awk '{printf "%+15.8E", \$2}')

    echo "dxs required: \$dx0 \$dy0 \$dz0"
    echo "dxs actually used: \$dx \$dy \$dz"

    # update to dxs actually used
    sed -i "s/dx(m).*/dx(m):              \$dx/"  $misfit_dir/diff_CMTSOLUTION
    sed -i "s/dy(m).*/dy(m):              \$dy/"  $misfit_dir/diff_CMTSOLUTION
    sed -i "s/dz(m).*/dz(m):              \$dz/"  $misfit_dir/diff_CMTSOLUTION

    ${python_exec} ${sem_utils_dir}/misfit/update_source_dxs.py \\
      "${db_file}" \\
      -- "\$dx" "\$dy" "\$dz"
    # -- means everything behind is positional,  so values like "-1.e-4" will be parsed correctly

  fi

done

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF

#====== search source parameters
cat <<EOF > $search_job
#!/bin/bash
#SBATCH -J ${event_id}.search
#SBATCH -o $search_job.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p $slurm_partition_cpu
#SBATCH -t $slurm_timelimit_search

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

cd $event_dir

# read derivatives of green's fuction
for tag in dxs dmt
do
  $python_exec $sem_utils_dir/misfit/read_perturbed_syn.py \\
    $db_file \\
    $event_dir/output_\${tag}/sac \\
    \${tag} \\
    --syn_is_grn
done

# grid search of source model
$python_exec $sem_utils_dir/misfit/grid_search_source.py \\
  $db_file \\
  $misfit_dir/grid_search_source.txt \\
  $misfit_dir/grid_search_source.pdf

# get optimal model
dxs_opt=\$(grep dxs_opt $misfit_dir/grid_search_source.txt | tail -n1 | awk '{print \$3}')
dmt_opt=\$(grep dmt_opt $misfit_dir/grid_search_source.txt | tail -n1 | awk '{print \$3}')
dt0_opt=\$(grep dt0_opt $misfit_dir/grid_search_source.txt | tail -n1 | awk '{print \$3}')
dtau_opt=\$(grep dtau_opt $misfit_dir/grid_search_source.txt | tail -n1 | awk '{print \$3}')

echo dxs_opt = \$dxs_opt
echo dmt_opt = \$dmt_opt
echo dt0_opt = \$dt0_opt
echo dtau_opt = \$dtau_opt

$python_exec $sem_utils_dir/misfit/make_updated_cmtsolution.py \\
  $db_file \\
  $misfit_dir/CMTSOLUTION.updated \\
  \$dt0_opt \$dtau_opt \$dxs_opt \$dmt_opt

echo
echo "Done: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

EOF
