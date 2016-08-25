#!/bin/bash

wkdir=$(pwd)

awk '$1!~/#/{print $2}' misfit.txt | sort -u |\
while read win
do
  echo $win
  output_fig=figure/${win}.pdf
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