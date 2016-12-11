#!/bin/bash

rm *job* *.list model mesh sem_utils

ln -sf ../model .
ln -sf ../mesh .
ln -sf ~/seiscode/sem_utils .
ln -sf sem_utils/utils/xsection/isc_d50km.txt .
ln -sf sem_utils/utils/xsection/fault_lines.txt .

# make list
sem_utils/utils/xsection/make_list_slice_gcircle.sh sem_utils slice_gcircle.list
sem_utils/utils/xsection/make_list_slice_sphere.sh slice_sphere.list

# make sbatch jobs
sem_utils/utils/xsection/make_sbatch_slice_gcircle.sh sem_utils mesh/DATABASES_MPI/ model slice_gcircle.list  dlnvs,kappa,eps,gamma nc slice_gcircle.job
sem_utils/utils/xsection/make_sbatch_slice_sphere.sh sem_utils mesh/DATABASES_MPI/ model slice_sphere.list  dlnvs,kappa,eps,gamma nc slice_sphere.job