#!/bin/bash

# copy files from previous iter_dir

source_dir=${1:?[arg]need previous iter_dir}

source_dir=$(readlink -f $source_dir)

cp -a $source_dir/*.txt .
cp -a $source_dir/*.job .
cp -a $source_dir/*.sh .
cp -a $source_dir/*.py .

find $source_dir -maxdepth 1 -type l | xargs -I@ cp -a @ .

rm model
ln -sf $source_dir/model_searched model