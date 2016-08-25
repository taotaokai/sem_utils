#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"
  cd $wkdir/$event_id/misfit

  for win in $(awk '$1!~/#/{print $2}' misfit.txt | sort -u )
  do
    echo $win
    output_fig=${win}.pdf
    input_figs=$(ls figure/*_${win}.pdf)
    err=$?
    if [ $err -ne 0 ]
    then
      echo "[ERROR] cannot find input figures"
    else
      rm $output_fig
      pdf_merge.sh $output_fig $input_figs
    fi
  done

done