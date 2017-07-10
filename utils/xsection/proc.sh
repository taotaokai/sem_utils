#!/bin/bash

wkdir=$(pwd)

sem_utils=~/seiscode/sem_utils

rm *job.o*

#rm *job* *.list model mesh sem_utils
#ln -sf ../model .
rm mesh; ln -sf ../mesh .
rm -rf model
mkdir model

ln -sf $(readlink -f ../model)/*_alpha.bin model/
ln -sf $(readlink -f ../model)/*_beta.bin model/
ln -sf $(readlink -f ../model)/*_phi.bin model/
ln -sf $(readlink -f ../model)/*_xi.bin model/
ln -sf $(readlink -f ../model)/*_eta.bin model/

ln -sf $(readlink -f ../mesh_REF/DATABASES_MPI)/*_reg1_vpv.bin model/
ln -sf $(readlink -f ../mesh_REF/DATABASES_MPI)/*_reg1_vsv.bin model/

cd model
ls *_vpv.bin | rename.sh "s/_vpv/_vp0/"
ls *_vsv.bin | rename.sh "s/_vsv/_vs0/"

ln -sf $sem_utils/utils/xsection/isc_d50km.txt .
ln -sf $sem_utils/utils/xsection/fault_lines.txt .

cd $wkdir
## make list
#sem_utils/utils/xsection/make_list_slice_gcircle.sh sem_utils slice_gcircle.list
#sem_utils/utils/xsection/make_list_slice_sphere.sh slice_sphere.list

# make sbatch jobs
$sem_utils/utils/xsection/make_sbatch_slice_gcircle.sh $sem_utils mesh/DATABASES_MPI/ model slice_gcircle.list  vp0,vs0,alpha,beta,phi,xi,eta nc slice_gcircle.job "ibrun"
$sem_utils/utils/xsection/make_sbatch_slice_sphere.sh $sem_utils mesh/DATABASES_MPI/ model slice_sphere.list vp0,vs0,alpha,beta,phi,xi,eta nc slice_sphere.job "ibrun"

#$sem_utils/utils/xsection/make_sbatch_slice_gcircle.sh $sem_utils mesh/DATABASES_MPI/ model slice_gcircle.list  alpha,beta,phi,xi,eta nc slice_gcircle.job "ibrun"
#$sem_utils/utils/xsection/make_sbatch_slice_sphere.sh $sem_utils mesh/DATABASES_MPI/ model slice_sphere.list alpha,beta,phi,xi,eta nc slice_sphere.job "ibrun"