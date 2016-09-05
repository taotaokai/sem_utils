#!/bin/bash

wkdir=$(pwd)

misfit_list=${1:?[arg]need misfit list (e.g. misfit.txt)}
figure_dir=${2:?[arg]need figure dir}
out_dir=${3:?[arg]need out_dir}

for win in $(awk '$1!~/#/{print $2}' $misfit_list | sort -u)
do
  echo ====== $win
  output_fig=$out_dir/${win}.pdf
  input_figs=$(ls $figure_dir/*_${win}.pdf)
  err=$?
  if [ $err -ne 0 ]
  then
    echo "[ERROR] cannot find input figures"
  else
    rm $output_fig
    pdf_merge.sh $output_fig $input_figs
  fi
done