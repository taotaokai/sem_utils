#!/bin/bash

# run a complete iteration

job_dep=${1:-1}

wkdir=$(pwd)

event_list=all_event.txt

mesh_dir=$wkdir/mesh
data_dir=$wkdir/events
utils_dir=$wkdir/utils
specfem_dir=$wkdir/specfem3d_globe
cmt_dir=$wkdir/CMTSOLUTION

base_dir=~/NEChina/compare_misfit_function/synthetic_inversion_with_1_event

rm *.job.o* *.job.submit *.log

## mesh
#ln -s $base_dir/iasp91/mesh .
##sbatch -d afterok:${job_dep} mesh_model.job | tee submit_jobs.log
#
## syn,misfit,kernel
##depjobs=$(grep "job_id=" submit_jobs.log | awk -F"=" '{printf "%s,",$2}' | sed "s/,$//")
#bash proc_event.sh $event_list syn,misfit,kernel ${job_dep} | tee submit_jobs.log

# post-proc
cd $wkdir
depjobs=$(grep "job_id=" submit_jobs.log | awk -F"=" '{printf "%s,",$2}' | sed "s/,$//")
bash post_proc_for_initial_iteration.sh $depjobs | tee post_proc_jobs.log

# perturb,search
cd $wkdir
depjobs=$(grep "job_id=" post_proc_jobs.log | awk -F"=" '{printf "%s,",$2}' | sed "s/,$//")
bash proc_event.sh $event_list perturb,search ${depjobs} | tee submit_jobs.log

# model_update
cd $wkdir
depjobs=$(grep "job_id=" submit_jobs.log | awk -F"=" '{printf "%s,",$2}' | sed "s/,$//")
sbatch -d afterok:${depjobs} model_update.job
