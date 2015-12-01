#!/bin/bash

# create new model from previous iterations
#   - model udpate direction: xsem_get_dmodel
#   - new model: xsem_add_dmodel

#====== command line args
control_file=${1:?must provide control_file}
mpi_exec=${2:-}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
fi
control_file=$(readlink -f $control_file)

# load parameters from control_file
source ${control_file}

echo
echo "Start updating model [$(date)]."
echo

#====== model update direction 
echo
echo "#====== get dmodel [$(date)]"
echo

if [ "${iter}" -le "${iter0}" ]
then

    echo "# no model update is needed [iter# $iter]"

elif [ "${iter_minus_one}" -eq "${iter0}" ] || [ "${use_lbfgs}" -eq 0 ]
then

    echo "use steepest descent [iter# $iter]"

    mkdir -p $model_dir/DATABASES_MPI
    cd $model_dir

    ${mpi_exec} \
        $sem_utils/bin/xsem_get_dmodel_steepest_descent \
        ${nproc} \
        ${prev_mesh_dir}/DATABASES_MPI \
        ${prev_kernel_dir}/DATABASES_MPI \
        "mu,lamda,rho" "_kernel" ${sd_scale_factor} 0 \
        ${model_dir}/DATABASES_MPI
else

    echo "use L-BFGS [iter# $iter]"

    mkdir -p $model_dir/DATABASES_MPI
    cd $model_dir

    # get dmodel,dkernel info from previous iteration steps
    n0=$(echo ${iter0} ${iter} ${nstep_lbfgs} |\
        awk '{n0=$2-$3; if(n0<=$1) n0=$1+1; print n0}')

    for i in $(seq ${iter_minus_one} -1 ${n0})
    do
        i=$(printf "%02d" $i)

        add_dmodel_log=$(readlink -f \
            ${base_dir}/iterations/iteration.${i}/model/xsem_add_dmodel.log)

        step_length=$(grep "^step_length" $add_dmodel_log | awk '{print $2}')

        dm_dir=$(readlink -f \
            ${base_dir}/iterations/iteration.${i}/model/DATABASES_MPI)

        dg_dir=$(readlink -f \
            ${base_dir}/iterations/iteration.${i}/kernel/DATABASES_MPI)

        printf "\"%s\"  \"%s\"  %g\n" ${dm_dir} ${dg_dir} $step_length

    done > ${model_dir}/dm_dg_alpha.list

    echo "# L-BFGS steps: ${n0}, ${iter_minus_one}"
    echo "# dmodel_dir dkernel_dir step_length:"
    echo
    cat ${model_dir}/dm_dg_alpha.list

    ${mpi_exec} \
       $sem_utils/bin/xsem_get_dmodel_lbfgs \
       ${nproc} \
       ${prev_mesh_dir}/DATABASES_MPI \
       ${prev_kernel_dir}/DATABASES_MPI \
       ${model_dir}/dm_dg_alpha.list \
       "mu,lamda,rho" "_kernel" 0 \
       ${model_dir}/DATABASES_MPI

fi


#====== create new model
echo
echo "#====== create new model [$(date)]"
echo

if [ "$iter" -le "${iter0}" ]
then
    mkdir -p $model_dir/DATABASES_MPI

    # link initial model
    cd $model_dir/DATABASES_MPI
    ln -sf $init_model_dir/DATABASES_MPI/proc*_reg1_v??.bin ./ 
    ln -sf $init_model_dir/DATABASES_MPI/proc*_reg1_eta.bin ./ 
    ln -sf $init_model_dir/DATABASES_MPI/proc*_reg1_rho.bin ./ 

elif [ -f ${prev_model_dir}/DATABASES_MPI/proc000000_reg1_vpv.bin ]
then
    mkdir -p $model_dir/DATABASES_MPI

    ${mpi_exec} $sem_utils/bin/xsem_add_dmodel_lamda_mu_to_tiso \
        $nproc \
        ${prev_mesh_dir}/DATABASES_MPI \
        ${prev_model_dir}/DATABASES_MPI \
        ${model_dir}/DATABASES_MPI \
        ${max_dlnv_allowed} \
        ${force_max_dlnv_allowed} \
        ${fix_rho} \
        ${model_dir}/DATABASES_MPI > ${model_dir}/xsem_add_dmodel.log
else
    echo "[ERROR] $prev_model_dir/DATABASES_MPI/proc000000_reg1_vpv.bin does NOT exist!"
    exit -1
fi

echo
echo "The model update is finished [$(date)]."
echo

#END
