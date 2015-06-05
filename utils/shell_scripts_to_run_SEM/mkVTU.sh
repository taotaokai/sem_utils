#!/bin/bash

#make VTU files of kernels
echo $(date)

# read Par_file
source Par_file

# region id
IREG=1

# each event directory for forward simulation
echo event directory ...
for evnm in $(cat $iter_dir/DATA/EVENTS)
do
    evnm_dir=$iter_dir/EVENTS/$evnm

    if [ ! -d $evnm_dir/VTU ];then
        mkdir $evnm_dir/VTU
    fi

    # combine volume data
    cd $evnm_dir
    seq 0 $((num_proc - 1)) > VTU/SLICES_ID
    topo_dir=$evnm_dir/DATABASES_MPI
    local_dir=$evnm_dir/DATABASES_MPI
    out_dir=$evnm_dir/VTU
    #for tag in alpha_kernel beta_kernel hess_kernel
    #for tag in beta_kernel alpha_kernel #hess_kernel
    #for tag in beta_kernel #alpha_kernel 
    for tag in mask_source #hess_kernel 
    do
        echo xcombine_vol_data
        $bin_dir/xcombine_vol_data VTU/SLICES_ID $tag $topo_dir $local_dir $out_dir 0 $IREG \
            2>&1 | tee $out_dir/$tag.log
    	echo combine volumen data successfully  
        # convert to vtu files
        echo mesh2vtu
        $bin_dir/mesh2vtu -i VTU/reg_${IREG}_${tag}.mesh -o VTU/reg${IREG}_${tag}.vtu \
            2>&1 | tee -a $out_dir/$tag.log
    	echo mesh2vtu successfully
    done
done

#END
