#!/bin/bash

#20150131 created
# setup folders for SEM simulation of one iteration 

# get parameters
source Par_file 

# check parameters 
if [ ! -d "$base_dir" ]
then
	echo "ERROR: base directory ($base_dir) does not exist!"
	exit
fi

if [ -z "$iter" ] || [ "$iter" -lt 0 ]
then
	echo "ERROR: iteration number ($iter) must equal or greater than 0"
	exit	
fi

if [ ! -f $config_dir/DATA/Par_file ]
then
	echo "ERROR: $config_dir/DATA/Par_file must exist!"
	exit
fi

if [ -d ${iter_dir} ]
then
	echo "WARN: iteration directory ($iter_dir) already exists!"
fi

if [ ! -d $init_model_dir ]
then
	echo "ERROR: initial model directory ($init_model_dir) does not exist"
	exit
fi

if [ $iter -gt 0 ] && [ ! -d $prev_iter_dir ]
then
	echo "ERROR: previous iteration directory ($prev_iter_dir) does not exist"
	exit
fi

#END
