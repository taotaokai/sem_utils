#!/bin/bash

wkdir=$(pwd)
utils=~/seiscode/sem_utils/utils/misfit_v0

event_list=${1:?[arg]need fdsnws-event list}
out_dir=${2:?[arg]need out dir}

chmod u+w -R $out_dir
rm -rf $out_dir
mkdir $out_dir

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  event_dir=$wkdir/$event_id

  cp $event_dir/misfit/misfit.txt $out_dir/${event_id}.txt

done

#cp $wkdir/misfit_par.py $out_dir/

cd $out_dir

# concantenate all misfit files 
tmp=$(mktemp -p .)
ls *.txt > list
for fname in $(cat list)
do
  cat $fname >> $tmp
done

$utils/plot_hist.py $tmp misfit_hist.pdf
rm $tmp