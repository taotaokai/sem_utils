#!/bin/bash

wkdir=$(pwd)
sem_utils=~/seiscode/sem_utils

mkdir grid_search_dmodel
ls */misfit/grid_search_dmodel.txt | awk -F"/" '{printf "cp %s grid_search_dmodel/%s.txt\n",$0,$1}' | bash

cd $wkdir/grid_search_dmodel
ls *.txt > list
$sem_utils/utils/misfit_v0/plot_grid_search_1d.py list
