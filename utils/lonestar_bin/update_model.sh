#!/bin/bash

# sum up all event kernels
#   - sum event kernels (cijkl, rho)
#   - reduce cijkl_kernel to (lamda,mu)_kernel
#   - get model udpate direction

#====== command line args
control_file=${1:?must provide control_file}
event_list=${2:?must provide event_id list}
mpi_exec=${3:-}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
fi
control_file=$(readlink -f $control_file)

if [ ! -f "$event_list" ]
then
    echo "[ERROR] invalid event_list: ", $event_list
fi
event_list=$(readlink -f $event_list)

# load parameters from control_file
source ${control_file}

echo "Start updating model [$(date)]."

#====== model gradient: sum up event kernels

#-- make event kernel mask
echo "#-- make event kernel mask [$(date)]"
for evid in $(grep -v ^# $event_list)
do
    echo $evid
    event_dir=$iter_dir/$evid
    # create source_xyz.list
    cd ${event_dir}
    src_vtk=$event_dir/OUTPUT_forward/source.vtk
    sed -n '/^POINTS/{n;p;}' $src_vtk > source_xyz.list
    cat source_xyz.list
    # create mask gll
    cd ${event_dir}
    ${mpi_exec} \
        $sem_utils/bin/xsem_make_source_depth_mask \
        $nproc ${mesh_dir}/DATABASES_MPI source_xyz.list \
        $source_gaussa $depth_pass $depth_gaussa \
        DATABASES_MPI
done

#-- make kernel_dir_list
cd ${iter_dir}
if [ ! -d model_update/DATABASES_MPI ]
then
    mkdir -p model_update/DATABASES_MPI 
fi

cd ${iter_dir}/model_update
awk '$1!~/^#/{printf "%s/%s/DATABASES_MPI\n", iter_dir, $1}' \
    iter_dir="${iter_dir}" $event_list > kernel_dir.list

cat kernel_dir.list

#-- sum up all event cijkl_kernel (use each event mask)
echo "#-- sum up all event cijkl_kernel [$(date)]"
cd ${iter_dir}/model_update
${mpi_exec} \
    $sem_utils/bin/xsem_sum_event_kernels_cijkl \
    $nproc $mesh_dir/DATABASES_MPI kernel_dir.list $use_mask $nroot_stack \
    DATABASES_MPI

#-- sum up all event rho_kernel (use each event mask)
echo "#-- sum up all event rho_kernel [$(date)]"
cd ${iter_dir}/model_update
${mpi_exec} \
    $sem_utils/bin/xsem_sum_event_kernels_1 \
    $nproc $mesh_dir/DATABASES_MPI kernel_dir.list "rho_kernel" \
    $use_mask $nroot_stack DATABASES_MPI

#-- reduce cijkl kernel to (lamda,mu)_kernel
echo "#-- reduce cijkl kernel to (lamda,mu)_kernel [$(date)]"
cd ${iter_dir}/model_update
${mpi_exec} \
    $sem_utils/bin/xsem_reduce_kernel_cijkl_to_lamda_mu \
    $nproc $mesh_dir/DATABASES_MPI DATABASES_MPI DATABASES_MPI

##-- kernel thresholding (lamda,mu,rho)_kernel
#echo "#-- kernel thresholding [$(date)]"
#cd ${iter_dir}/model_update
#for ker_name in lamda_kernel mu_kernel rho_kernel
#do
#    ${mpi_exec} $sem_utils/bin/xsem_pdf \
#        $nproc $mesh_dir/DATABASES_MPI DATABASES_MPI $ker_name  \
#        1000 1 ${ker_name}_abs_pdf.txt
#
#    zc=$(awk '$1!~/#/{a+=$3; if(a>cutoff){print $2; exit}}' \
#        cutoff=$threshold_corner ${ker_name}_abs_pdf.txt)
#
#    echo "### ${ker_name}: corner amplitude is $zc"
#
#    ${mpi_exec} $sem_utils/bin/xsem_thresholding \
#        $nproc $mesh_dir/DATABASES_MPI DATABASES_MPI $ker_name  \
#        $zc $threshold_rmax DATABASES_MPI ${ker_name}_precond
#done

#====== model update direction 

##-- create source_xyz.list
#cd ${iter_dir}/model_update
#for evid in $(grep -v ^# $event_list)
#do
#    event_dir=$iter_dir/$evid
#    src_vtk=$event_dir/OUTPUT_forward/source.vtk
#    sed -n '/^POINTS/{n;p;}' $src_vtk
#done > source_xyz.list
# 
##-- create mask gll
#cd ${iter_dir}/model_update
#${mpi_exec} \
#    $sem_utils/bin/xsem_make_source_depth_mask \
#    $nproc $mesh_dir/DATABASES_MPI source_xyz.list \
#    $source_mask_radius $stop_depth $pass_depth \
#    DATABASES_MPI

#-- get dkernel
echo "#-- get dkernel [$(date)]"
cd ${iter_dir}/model_update
if [ "${iter}" -gt 1 ]
then
    echo "get model gradient update"
    $sem_utils/bin/xsem_math \
        $nproc $mesh_dir/DATABASES_MPI \
        DATABASES_MPI "mu_kernel,lamda_kernel,rho_kernel" \
        $prev_iter_dir/model_update/DATABASES_MPI "mu_kernel,lamda_kernel,rho_kernel" \
        sub \
        DATABASES_MPI "mu_dkernel,lamda_dkernel,rho_dkernel"
else
    echo "first iteration, no dkernel"
fi

#-- get model update direction
echo "#-- get dmodel [$(date)]"
cd ${iter_dir}/model_update
if [ "${iter}" -le 1 ] || [ "${use_lbfgs}" -eq 0 ]
then 
    echo "use steepest descent"

    ${mpi_exec} \
        $sem_utils/bin/xsem_get_dmodel_steepest_descent \
        $nproc $mesh_dir/DATABASES_MPI DATABASES_MPI \
        "mu,lamda,rho" "_kernel" ${sd_scale_factor} 0 \
        DATABASES_MPI
else
    echo "use l-BFGS"

    # get dmodel,dkernel info from previous iteration steps
    n0=$(echo ${iter} ${nstep_lbfgs} | \
        awk '{n0=$1-$2; if(n0<2) n0=2; print n0}')
    for i in $(seq ${iter} -1 ${n0})
    do
        add_dmodel_log=$(printf \
            "$wkdir/iterations/iteration.%02d/mesh/xsem_add_dmodel.log" $i)
        step_length=$(grep "step_length:" $add_dmodel_log | awk '{print $2}')
        printf "\"%s/iterations/iteration.%02d/model_update/DATABASES_MPI\" \
\"%s/iterations/iteration.%02d/model_update/DATABASES_MPI\" \
%f\n" $wkdir $((i-1)) $wkdir $i $step_length
    done > dm_dg_alpha.list

   ${mpi_exec} \
       $sem_utils/bin/xsem_get_dmodel_lbfgs \
       $nproc $mesh_dir/DATABASES_MPI DATABASES_MPI \
       dm_dg_alpha.list "mu,lamda,rho" "_kernel" 0 \
       DATABASES_MPI
fi

echo "The model update is finished [$(date)]."
