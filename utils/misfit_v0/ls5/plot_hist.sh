#!/bin/bash

wkdir=$(pwd)
utils_dir=/home1/03244/ktao/seiscode/sem_utils/utils/misfit_v0

file_list=${1:?[arg]need misfit file list}
pdf_file=${2:?[arg]need output pdf name}

#====== concantenate all misfit files 
tmp=$(mktemp -p .)

for fname in $(cat $file_list)
do
  cat $fname >> $tmp
done

sed -i 's/Rayleigh/surface/;s/Love/surface/' $tmp

#====== make histograms of dt and cc
$utils_dir/hist_misfit.py $tmp $pdf_file

#clear
rm $tmp