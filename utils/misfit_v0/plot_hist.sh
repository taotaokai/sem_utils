#!/bin/bash

wkdir=$(pwd)
utils_dir=/home1/03244/ktao/seiscode/sem_utils/utils/misfit_v0

event_list=${1:?[arg]need event_list}
pdf_file=${2:?[arg]need output pdf name}

#====== concantenate all misfit files 
tmp=$(mktemp -p .)
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do
  misfit_file=$wkdir/$event_id/misfit/misfit.txt
  cat $misfit_file >> $tmp
done
sed -i 's/Rayleigh/surface/;s/Love/surface/' $tmp

#====== make histograms of dt and cc
$utils_dir/hist_misfit.py $tmp $pdf_file

#clear
rm $tmp