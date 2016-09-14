#!/bin/bash

log_file=${1:?[arg]need .o12345 file}

for model in dlnvs kappa eps gamma
do

  minval=$(grep $model $log_file | awk '{printf "%f\n",$4}' | sort -nr | tail -n1)
  maxval=$(grep $model $log_file | awk '{printf "%f\n",$5}' | sort -n  | tail -n1)
  echo $model $minval $maxval

done