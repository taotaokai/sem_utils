#!/bin/bash

rm *job* *.list model mesh sem_utils

ln -sf ../model .
ln -sf ../mesh .
ln -sf ~/seiscode/sem_utils .

# make list
sem_utils/utils/xsection/make_list_slice_gcircle.sh sem_utils slice_gcircle.list
sem_utils/utils/xsection/make_list_slice_sphere.sh slice_sphere.list

# make sbatch jobs
sem_utils/utils/xsection/make_sbatch_slice_gcircle.sh sem_utils mesh/DATABASES_MPI/ model slice_gcircle.list  vsv,vsh,vpv,vph,rho,eta nc slice_gcircle.job
sem_utils/utils/xsection/make_sbatch_slice_sphere.sh sem_utils mesh/DATABASES_MPI/ model slice_sphere.list  vsv,vsh,vpv,vph,rho,eta nc slice_sphere.job