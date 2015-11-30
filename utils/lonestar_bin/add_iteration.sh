#!/bin/bash

# add new iteration, create corresponding directories:
# iteration.xx/

#====== command line args
control_file=${1:?must provide control_file}

# source parameters in control_file
source ${control_file}

#====== create iteration directory
if [ ! -d ${iter_dir} ];then
    mkdir -p $iter_dir
else
    echo "[WARNING] iter_dir=$iter_dir already exits!"
    exit -1
fi

cd $wkdir
if [ ! -d iterations ]
then
    mkdir iterations
fi
cd iterations
ln -sf $iter_dir
