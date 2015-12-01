#!/bin/bash

# submit job scripts for one iteration:
#   - model update
#   - mesh
#   - forward+adjoint for each event 
#   - model gradient, preconditioner, direction 

control_file=${1:?must provide control_file}
event_list=${2:?must provide event_id_list}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
    exit -1
fi
control_file=$(readlink -f $control_file)

if [ ! -f "$event_list" ]
then
    echo "[ERROR] invalid event_id list: ", $event_list
    exit -1
fi
event_list=$(readlink -f $event_list)

# load parameters in control_file
source ${control_file}

#====== iteration
echo
echo "====== iteration: ${iter}"
echo

if [ ! -d ${iter_dir} ];then
    mkdir -p $iter_dir
else
    echo "[WARNING] iter_dir=$iter_dir already exits!"
    exit -1
fi

cd $base_dir
if [ ! -d iterations ]
then
    mkdir iterations
fi
cd iterations
ln -sf $iter_dir

#====== mesh
echo
echo ====== mesh: submit job script
echo

cd $base_dir
job_file=${iter_dir}/mesh.job
cat <<EOF > ${job_file}
#!/bin/bash
#$ -V                                       # Inherit the submission environment 
#$ -cwd                                     # Start job in submission directory
#$ -N mesh.$iter                            # Job Name
#$ -j y                                     # combine stderr & stdout into stdout  
#$ -o ${iter_dir}/mesh.${iter}.o\$JOB_ID    # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $nproc_request                 # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                                # Queue name
#$ -l h_rt=$run_time_mesher                 # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu                    # email 
#$ -m bea                                   # email info: begin/end/abort
#$ -hold_jid -1                             # dependent job id

export MY_NSLOTS=$nproc

echo
echo "Job begins: JOB_ID=\${JOB_ID} [\$(date)]"
echo

echo
echo "====== create new model [\$(date)]"
echo
$base_dir/bin/update_model.sh $control_file ibrun

echo
echo "====== mesh [\$(date)]"
echo
echo "# setup mesh directory"
$base_dir/bin/setup_mesh.sh $control_file
echo
echo "# run mesher"
cd $mesh_dir
ibrun $build_dir/bin/xmeshfem3D

echo
echo "Job ends: JOB_ID=\${JOB_ID} [\$(date)]"
echo

EOF

qsub ${job_file} > ${job_file}_qsub
if [ "$?" -ne 0 ]
then
    echo "[ERROR] failed to qsub ${job_file}"
    exit -1
fi

# get mesh job_id
mesh_jid=$(awk 'END{print $3}' ${job_file}_qsub)
echo "mesh job id is $mesh_jid"

#====== forward+adjoint
echo
echo ====== events: submit job script
echo

cd $base_dir
event_jid_list=""
for evid in $(grep -v ^# $event_list)
do
    event_dir=${iter_dir}/$evid
    job_file=${iter_dir}/${evid}.job

cat<<EOF > $job_file
#!/bin/bash
#$ -V                                       # Inherit the submission environment 
#$ -cwd                                     # Start job in submission directory
#$ -N $evid.$iter                          # Job Name
#$ -j y                                     # combine stderr & stdout into stdout  
#$ -o ${iter_dir}/${evid}.${iter}.o\$JOB_ID # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $nproc_request                 # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                                # Queue name
#$ -l h_rt=$run_time_event                  # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu                    # email 
#$ -m bea                                   # email info: begin/end/abort
#$ -hold_jid ${mesh_jid:--1}                # dependent job id

export MY_NSLOTS=$nproc

echo
echo "Job begins: JOB_ID=\${JOB_ID} [\$(date)]"
echo

echo
echo "====== forward [\$(date)]"
echo
$base_dir/bin/setup_event.sh $control_file $evid
cd $iter_dir/$evid
ibrun $build_dir/bin/xspecfem3D

echo
echo "====== measure misfit [\$(date)]"
echo
$base_dir/bin/measure_adjoint.sh $control_file $evid

echo
echo "====== adjoint(model) [\$(date)]"
echo
$base_dir/bin/setup_adjoint.sh $control_file $evid
cd $iter_dir/$evid
ibrun $build_dir/bin/xspecfem3D

echo
echo "Job ends: JOB_ID=\${JOB_ID} [\$(date)]"
echo

EOF

    qsub $job_file > ${job_file}_qsub 
    if [ "$?" -ne 0 ]
    then
        echo "[ERROR] failed to qsub ${job_file}"
        exit -1
    fi 

    event_jid=$(awk 'END{print $3}' ${job_file}_qsub)
    if [ -z "${event_jid_list}" ]
    then
        event_jid_list=${event_jid}
    else
        event_jid_list=${event_jid_list},${event_jid}
    fi
done

echo "event job id's: $event_jid_list"

#====== get model gradient 
echo
echo ====== update: submit job script
echo

cd $base_dir
job_file=${iter_dir}/kernel.job

cat <<EOF > ${job_file}
#!/bin/bash
#$ -V                                       # Inherit the submission environment 
#$ -cwd                                     # Start job in submission directory
#$ -N update.$iter                          # Job Name
#$ -j y                                     # combine stderr & stdout into stdout  
#$ -o ${iter_dir}/kernel.${iter}.o\$JOB_ID  # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12                             # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                                # Queue name
#$ -l h_rt=$run_time_update                 # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu                    # email 
#$ -m bea                                   # email info: begin/end/abort
#$ -hold_jid $event_jid_list                # dependent job id

echo
echo "Job begins: JOB_ID=\${JOB_ID} [\$(date)]"
echo

echo
echo "====== update model gradient [\$(date)]"
echo

$base_dir/bin/update_kernel.sh $control_file $event_list ibrun

echo
echo "Job ends: JOB_ID=\${JOB_ID} [\$(date)]"
echo

EOF

qsub ${job_file} > ${job_file}_qsub
if [ "$?" -ne 0 ]
then
    echo "[ERROR] failed to qsub ${job_file}"
    exit -1
fi

#END
