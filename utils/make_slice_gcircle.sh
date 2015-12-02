#!/bin/bash

# interpolate SEM gll models to create gcircle slices

#====== command line args
control_file=${1:?[arg] need control_file}
event_id=${2:?[arg] need evid}
slice_list=${3:?[arg] need slice_list}
mpi_exec=${4:?[arg] need mpi_exec}
model_names=${5:-vsv,vsh,vpv,vph,rho,eta}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: "  $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)
# load parameters in control_file
source ${control_file}

if [ ! -f "$slice_list" ]
then
    echo "[ERROR] invalid slice_list: "  $slice_list
    exit -1
fi
slice_list=$(readlink -f $slice_list)

event_dir=${iter_dir}/${event_id}
if [ ! -d "${event_dir}" ]
then
    echo "[ERROR] <event_id> directory deos not exist: "  $event_id
    exit -1
fi

#====== interpolate SEM gll models

mkdir -p $event_dir/xsection

#-- create job script
job_file=${event_dir}/slice_gcircle.job

cat<<EOF > $job_file
#!/bin/bash
#$ -V                              # Inherit the submission environment 
#$ -cwd                            # Start job in submission directory
#$ -N slice_gcircle.${iter}             # Job Name
#$ -j y                            # combine stderr & stdout into stdout  
#$ -o ${event_dir}/slice_gcircle.${iter}.o\$JOB_ID    # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12                    # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                       # Queue name
#$ -l h_rt=00:30:00                # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu           # email 
#$ -m bea                          # email info: begin/end/abort
#$ -hold_jid -1                    # dependent job id

EOF

grep -v "^#" $slice_list |\
while read lat0 lon0 azimuth theta0 theta1 ntheta r0 r1 nr fname
do
    out_file=${fname}.nc

cat<<EOF >> $job_file
echo
echo \$(date)
echo "# $fname: $lat0 $lon0 $azimuth $theta0 $theta1 $ntheta $r0 $r1 $nr "
echo
${mpi_exec} \
    $sem_utils/bin/xsem_slice_gcircle \
    $nproc \
    $mesh_dir/DATABASES_MPI \
    $event_dir/DATABASES_MPI \
    $model_names \
    $lat0 $lon0 $azimuth \
    $theta0 $theta1 $ntheta \
    $r0 $r1 $nr \
    $event_dir/xsection/$out_file

EOF

done
