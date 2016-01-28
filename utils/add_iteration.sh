#!/bin/bash

# add new iteration, create corresponding directories:
# iteration.xx/

#====== command line args
control_file=${1:?[arg] need control_file}

# check inputs
if [ ! -f "$control_file" ]
then
    echo "[ERROR] invalid control_file: ", $control_file
    exit 1
fi
control_file=$(readlink -f $control_file)

# load parameters in control_file
source ${control_file}

#====== create iteration directory
if [ ! -d "${iter_dir}" ]
then
    mkdir -p $iter_dir
else
    echo "[ERROR] iter_dir=$iter_dir already exits!"
    exit 1
fi

cd $base_dir
if [ ! -d iterations ]
then
    mkdir iterations
fi

cd iterations
ln -sf $iter_dir
