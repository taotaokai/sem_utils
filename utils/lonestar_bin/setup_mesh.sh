#!/bin/bash

# setup mesh folders, generate the batch script to run SEM meshfem3D


#====== command line args
control_file=${1:?must provide control_file}
mpi_exec=${2:-}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
fi
control_file=$(readlink -f $control_file)

# source parameters in control_file
source ${control_file}

#====== create mesh dir
if [ ! -d ${mesh_dir} ];then
    mkdir -p $mesh_dir
else
    echo "[WARNING] mesh_dir=$mesh_dir already exits!"
    exit -1
fi

#====== setup mesh dir
cd $mesh_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

# create necessary files in mesh/DATA/
cd $mesh_dir/DATA
# all data files: topography, bathymetry, etc.
ln -sf $build_dir/DATA/* ./
# modify Par_file
rm Par_file
cp -L $config_dir/DATA/Par_file .
if [ ${iter} -eq 0 ] # save mesh file for starting model
then
    sed -i "/^SAVE_MESH_FILES/s/=.*/= .true./" Par_file
else
    sed -i "/^SAVE_MESH_FILES/s/=.*/= .false./" Par_file
fi

# set model
if [ "$iter" -eq 0 ]
then
    cd $mesh_dir/DATA
    ln -sf $config_dir/starting_model GLL
elif [ -d ${prev_iter_dir}/model_update/DATABASES_MPI ]
then
    # get new model
    cd $mesh_dir/DATABASES_MPI
    ln -sf ${prev_iter_dir}/model_update/DATABASES_MPI/*_dmodel.bin ./
    ${mpi_exec} $sem_utils/bin/xsem_add_dmodel_lamda_mu_to_tiso \
        $nproc \
        ${prev_iter_dir}/mesh/DATABASES_MPI \
        ${prev_iter_dir}/mesh/DATABASES_MPI \
        ${prev_iter_dir}/model_update/DATABASES_MPI \
        ${max_dlnv_allowed} \
        ${force_max_dlnv_allowed} \
        ${fix_rho} \
        ${mesh_dir}/DATABASES_MPI > ${mesh_dir}/xsem_add_dmodel.log
    # link model dir
    cd $mesh_dir/DATA
    ln -sf $mesh_dir/DATABASES_MPI GLL
else
    echo "[ERROR] $prev_iter_dir/model_update/DATABASES_MPI does NOT exist!"
    exit -1
fi

# backup parameter files into OUTPUT_FILES
cd $mesh_dir
cp -L DATA/Par_file OUTPUT_FILES
cp -L DATA/CMTSOLUTION OUTPUT_FILES
cp -L DATA/STATIONS OUTPUT_FILES

#====== make job script
#cat <<EOF > $mesh_dir/mesher.job
##$ -V                              # Inherit the submission environment 
##$ -cwd                            # Start job in submission directory
##$ -N mesher.$iter                 # Job Name
##$ -j y                            # combine stderr & stdout into stdout  
##$ -o $mesh_dir/mesher.o           # Name of the output file (eg. myMPI.oJobID)
##$ -pe 12way $nproc_request        # Requests 12 cores/node, 24 cores total: 12way 24
##$ -q normal                       # Queue name
##$ -l h_rt=$run_time_mesher        # Run time (hh:mm:ss) - 1.5 hours
##$ -M kai.tao@utexas.edu           # email 
##$ -m bea                          # email info: begin/end/abort
##$ -hold_jid -1                    # dependent job id
#
#export MY_NSLOTS=$nproc
#
#cd $mesh_dir
#ibrun $build_dir/bin/xmeshfem3D
#
#EOF
#
#echo "The mesher can be run as"
#echo ">> qsub $mesh_dir/mesher.job"
