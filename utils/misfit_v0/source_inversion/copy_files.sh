#!/bin/bash

# copy files from previous iter_dir

source_dir=${1:?[arg]need previous iter_dir}

source_dir=$(readlink -f $source_dir)

find $source_dir -maxdepth 1 -type l | xargs -I@ cp -a @ .

rm CMTSOLUTION_initial
ln -s $(readlink -f $source_dir/CMTSOLUTION_updated) CMTSOLUTION_initial
#ln -s $(readlink -f $source_dir/mesh) mesh

rm mesh; ln -s $(readlink -f $source_dir/mesh) mesh
cp -a $source_dir/proc.sh .
cp -a $source_dir/to_run.txt .
#cp -a $source_dir/misfit_par.py .
cp -a $source_dir/misfit_par .

