#!/bin/bash

# setup event folder, generate the batch script to run SEM xspecfem3D
# for forward simulation

#====== command line args
control_file=${1:-control_file}

# source parameters in control_file
if [ ! -f $control_file ]
then
    echo "[ERROR] control file doesn't exit! [$control_file]"
    exit -1
fi
source ${control_file}

#====== create event dir
cd $mesh_dir
if [ ! -d cross_sections ]
then
    mkdir cross_sections
fi
cd cross_sections
ln -sf $config_dir/cross_sections/*.xyz .

#====== make job script
job_file=$mesh_dir/interp.job
cat <<EOF > $job_file
#$ -V                              # Inherit the submission environment 
#$ -cwd                            # Start job in submission directory
#$ -N $evid.i$iter                 # Job Name
#$ -j y                            # combine stderr & stdout into stdout  
#$ -o $mesh_dir/interp.o          # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12                    # Requests 12 cores/node, 24 cores total: 12way 24
#$ -q normal                       # Queue name
#$ -l h_rt=$run_time_interp        # Run time (hh:mm:ss) - 1.5 hours
#$ -M kai.tao@utexas.edu           # email 
#$ -m bea                          # email info: begin/end/abort
#$ -hold_jid -1                    # dependent job id

EOF

cd $mesh_dir/cross_sections

for xyz_file in $(ls *.xyz)
do
    out_file=${xyz_file%.xyz}.out

cat<<EOF >> $job_file

cd $mesh_dir
$sem_utils/bin/xsem_interp_xyz \
     DATABASES_MPI $nproc DATA/GLL vpv,vph,vsv,vsh,eta,rho \
     cross_sections/$xyz_file cross_sections/$out_file

EOF

done

#
echo "The interpolation can be run as"
echo "qsub $mesh_dir/interp.job"
