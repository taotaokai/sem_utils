#!/bin/bash

# copy files from previous iter_dir

source_dir=${1:?[arg]need previous iter_dir}
target_dir=${2:?[arg]need current iter_dir}

source_dir=$(readlink -f $source_dir)

mkdir $target_dir

#cp -a $source_dir/old .
cp -a $source_dir/*.txt $target_dir/
cp -a $source_dir/*.job $target_dir/
cp -a $source_dir/*.sh $target_dir/
cp -a $source_dir/*.py $target_dir/
cp -a $source_dir/misfit_par $target_dir/
cp -a $source_dir/control_file $target_dir/

find $source_dir -maxdepth 1 -type l | xargs -I@ cp -a @ $target_dir/

rm $target_dir/model
ln -sf $source_dir/model_updated $target_dir/model