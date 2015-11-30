#!/bin/bash

# submit job scripts for one iteration:
#   - model update
#   - mesh
#   - forward+adjoint for each event 
#   - model gradient, preconditioner, direction 

control_file=${1:?must provide control_file}
event_list=${2:?must provide event_id_list}

control_file=$(readlink -f $control_file)
source ${control_file}

#====== iteration 
echo
echo "====== iteration: ${iter}"
echo
$wkdir/bin/add_iteration.sh $control_file

#====== mesh
echo
echo ====== mesh: submit job script
echo
cat <<EOF > ${iter_dir}/mesher.job
#$ -V                                 # Inherit the submission environment 
#$ -cwd                               # Start job in submission directory
#$ -N mesher.$iter                    # Job Name
#$ -j y                               # combine stderr & stdout into stdout  
#$ -o ${iter_dir}/mesher.\${JOB_ID}.o # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $nproc_request           # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                          # Queue name
#$ -l h_rt=$run_time_mesher           # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu              # email 
#$ -m bea                             # email info: begin/end/abort
#$ -hold_jid -1                       # dependent job id

export MY_NSLOTS=$nproc

echo
echo "Job begins: JOB_ID=\${JOB_ID} [\$(date)]"
echo

echo
echo "====== mesh [\$(date)]"
echo
$wkdir/bin/setup_mesh.sh $control_file
cd $mesh_dir
ibrun $build_dir/bin/xmeshfem3D

echo
echo "Job ends: JOB_ID=\${JOB_ID} [\$(date)]"
echo

EOF

qsub ${iter_dir}/mesher.job > ${iter_dir}/mesher.qsub

# get mesh job_id
mesh_jid=$(awk 'END{print $3}' ${iter_dir}/mesher.qsub)
echo "mesh job id is $mesh_jid"

#====== forward+adjoint
echo
echo ====== events: submit job script
echo

event_jid_list=""
for evid in $(grep -v ^# $event_list)
do
    event_dir=${iter_dir}/$evid
    job_file=${iter_dir}/${evid}.job

cat<<EOF > $job_file
#$ -V                                  # Inherit the submission environment 
#$ -cwd                                # Start job in submission directory
#$ -N $evid.f$iter                     # Job Name
#$ -j y                                # combine stderr & stdout into stdout  
#$ -o ${iter_dir}/${evid}.\${JOB_ID}.o # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $nproc_request            # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                           # Queue name
#$ -l h_rt=$run_time_event             # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu               # email 
#$ -m bea                              # email info: begin/end/abort
#$ -hold_jid ${mesh_jid:--1}           # dependent job id

export MY_NSLOTS=$nproc

echo
echo "Job begins: JOB_ID=\${JOB_ID} [\$(date)]"
echo

echo
echo "====== forward [\$(date)]"
echo
$wkdir/bin/setup_event.sh $control_file $evid
cd $iter_dir/$evid
ibrun $build_dir/bin/xspecfem3D

echo
echo "====== measure misfit [\$(date)]"
echo
$wkdir/bin/measure_adjoint.sh $control_file $evid

echo
echo "====== adjoint(model) [\$(date)]"
echo
$wkdir/bin/setup_adjoint.sh $control_file $evid
cd $iter_dir/$evid
ibrun $build_dir/bin/xspecfem3D

echo
echo "Job ends: JOB_ID=\${JOB_ID} [\$(date)]"
echo

EOF

    qsub $job_file > ${iter_dir}/${evid}.qsub 
    event_jid=$(awk 'END{print $3}' ${iter_dir}/${evid}.qsub)
    if [ -z "${event_jid_list}" ]
    then
        event_jid_list=${event_jid}
    else
        event_jid_list=${event_jid_list},${event_jid}
    fi
done

echo "event job id's: $event_jid_list"

#====== get model update direction
echo
echo ====== update: submit job script
echo
cat <<EOF > ${iter_dir}/update.job
#$ -V                                 # Inherit the submission environment 
#$ -cwd                               # Start job in submission directory
#$ -N update.$iter                    # Job Name
#$ -j y                               # combine stderr & stdout into stdout  
#$ -o ${iter_dir}/update.\${JOB_ID}.o # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $nproc_request           # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                          # Queue name
#$ -l h_rt=$run_time_update           # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu              # email 
#$ -m bea                             # email info: begin/end/abort
#$ -hold_jid $event_jid_list          # dependent job id

export MY_NSLOTS=$nproc

echo
echo "Job begins: JOB_ID=\${JOB_ID} [\$(date)]"
echo

echo
echo "====== update [\$(date)]"
echo

$wkdir/bin/update_model.sh $control_file $event_list ibrun

echo
echo "Job ends: JOB_ID=\${JOB_ID} [\$(date)]"
echo

EOF

qsub ${iter_dir}/update.job > ${iter_dir}/update.qsub


#echo "Job finised. [\$(date)]"
