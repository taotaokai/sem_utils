#!/bin/bash

wkdir=$(pwd)
sem_utils=/home1/03244/ktao/seiscode/sem_utils

#mkdir grid_search_dvsh
#ls */misfit/grid_search_dvsh.txt | awk -F"/" '{printf "cp %s grid_search_dvsh/%s.txt\n",$0,$1}' | bash
#cd $wkdir/grid_search_dvsh
#ls *.txt > list
#$sem_utils/utils/misfit_v0/ls5/plot_grid_search_1d.py list
#
#
#cd $wkdir
#mkdir grid_search_dvp_dvsv
#ls */misfit/grid_search_dvp_dvsv.txt | awk -F"/" '{printf "cp %s grid_search_dvp_dvsv/%s.txt\n",$0,$1}' | bash
#cp misfit_par.py grid_search_dvp_dvsv/
#cd $wkdir/grid_search_dvp_dvsv
#ls *.txt > list
#$sem_utils/utils/misfit_v0/ls5/plot_grid_search_dvp_dvsv.py misfit_par.py list


mkdir grid_search_dmodel
ls */misfit/grid_search_dmodel.txt | awk -F"/" '{printf "cp %s grid_search_dmodel/%s.txt\n",$0,$1}' | bash
#cp misfit_par.py grid_search_dmodel/

cd $wkdir/grid_search_dmodel

ls *.txt > list
$sem_utils/utils/misfit_v0/ls5/plot_grid_search_1d.py list
