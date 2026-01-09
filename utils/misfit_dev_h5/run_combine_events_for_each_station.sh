#!/bin/bash

awk 'NR!=1' misfit_surf_0.05Hz.csv | sort -t, -k5 -k8 -k11 -n |\
while read -r line
do
  net=$(echo $line | awk -F, '{print $1}')
  sta=$(echo $line | awk -F, '{print $2}')
  echo "echo ===== $net $sta \($line\)"
  evdir=stations/${net}_${sta}
  h5file=$evdir/${net}_${sta}.h5
  echo mkdir -p $evdir
  echo python sem_utils/misfit/combine_misfit_for_one_station.py misfit_h5.lst $net $sta ${h5file}
  for win in Surf_Z_40-80sec Surf_T_40-80sec Surf_R_40-80sec Surf_Z_20-40sec Surf_T_20-40sec Surf_R_20-40sec Pwin_Z_10-100sec
  do
    fig=$evdir/${net}_${sta}_${win}.pdf
    echo python sem_utils/misfit/plot_misfit_for_one_station.py ${h5file} ${win} ${fig}
  done
done
