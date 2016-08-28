#!/bin/bash

mkdir misfit

ls */misfit/misfit.txt | awk -F"/" '{printf "cp %s misfit/%s.txt\n",$0,$1}' | bash

cp utils/ls5/plot_hist.sh misfit/

cd misfit
ls *.txt > list
./plot_hist.sh list misfit_hist.pdf