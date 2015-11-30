#!/bin/bash

control_file=${1:-control_file}
evid=C072703C

control_file=$(readlink -f $control_file)
source ${control_file}

#====== create job file
cat<<EOF
#$ -V                         # Inherit the submission environment 
#$ -cwd                       # Start job in submission directory
#$ -N doall.$iter             # Job Name
#$ -j y                       # combine stderr & stdout into stdout  
#$ -o doall.$iter.o           # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way $nproc_request   # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                  # Queue name
#$ -l h_rt=$run_time_all      # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu      # email 
#$ -m bea                     # email info: begin/end/abort
#$ -hold_jid -1               # dependent job id

export MY_NSLOTS=$nproc

echo "Job started. [\$(date)]"

$wkdir/bin/add_iteration.sh $control_file

$wkdir/bin/setup_mesh.sh $control_file
cd $mesh_dir
ibrun $build_dir/bin/xmeshfem3D

$wkdir/bin/setup_event.sh $control_file $evid
cd $iter_dir/$evid
ibrun $build_dir/bin/xspecfem3D

$wkdir/bin/measure_adjoint.sh $control_file $evid
$wkdir/bin/setup_adjoint.sh $control_file $evid
cd $iter_dir/$evid
ibrun $build_dir/bin/xspecfem3D

$wkdir/bin/update_model_1event.sh $control_file $evid

echo "Job finised. [\$(date)]"

EOF
