#!/bin/bash

wkdir=$(pwd)

prev_dir=${1:?[arg]need previous iteration dir}
curr_dir=${2:?[arg]need current iteration dir}
job_dep=${3:--1}

mkdir -p $wkdir/$curr_dir

cd $wkdir/$curr_dir
$wkdir/$prev_dir/copy_files.sh $wkdir/$prev_dir/
