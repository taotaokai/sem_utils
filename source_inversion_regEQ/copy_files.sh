#!/bin/bash

# copy files from previous iter_dir

prev_dir=${1:?[arg]need previous iter_dir}
current_dir=${2:?[arg]need current iter_dir}

prev_dir=$(readlink -f $prev_dir)
current_dir=$(readlink -f $current_dir)

if [ x$prev_dir == x$current_dir ]
then
  echo "[ERROR] prev_dir and current_dir cannot be the same directory!" 
  exit -1
fi

if [ ! -d "$current_dir" ]
then
  echo "[WARN] create current_dir ($current_dir)"
  mkdir -p $current_dir
fi

find $prev_dir -maxdepth 1 -type l | xargs -I@ cp -a @ $current_dir

cd $current_dir

rm CMTSOLUTION_initial
#ln -s $(readlink -f $prev_dir/CMTSOLUTION_updated) CMTSOLUTION_initial
cp -rL $(readlink -f $prev_dir/CMTSOLUTION_updated) CMTSOLUTION_initial
#ln -s $(readlink -f $prev_dir/mesh) mesh

rm mesh; ln -s $(readlink -f $prev_dir/mesh) mesh
#cp -a $prev_dir/proc.sh .
#cp -a $prev_dir/to_run.txt .
#cp -a $prev_dir/event.txt .
cp -a $prev_dir/*.txt .
#cp -a $prev_dir/misfit_par.py .
cp -a $prev_dir/misfit_par .
cp -a $prev_dir/control_file* .

echo !!!! modify control_file
