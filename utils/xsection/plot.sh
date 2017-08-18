#!/bin/bash
#SBATCH -J plot
#SBATCH -o plot.job.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 01:30:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

sem_utils=~/seiscode/sem_utils

#rm -rf figure/*

title=stage10.iter04

mkdir figure_alpha_beta figure_phi_xi_eta
awk '$1!~/#/&&$7<=220' slice_sphere.list > slice_sphere_220km.list

srun -n1 $sem_utils/utils/xsection/plot_slice_gcircle_alpha_beta_phi_xi_eta.sh nc slice_gcircle.list ${title} xi,beta,alpha figure_alpha_beta &
srun -n1 $sem_utils/utils/xsection/plot_slice_sphere_alpha_beta_phi_xi_eta.sh nc slice_sphere.list ${title} beta,alpha figure_alpha_beta &
#srun -n1 $sem_utils/utils/xsection/plot_slice_gcircle_alpha_beta_phi_xi_eta.sh nc slice_gcircle.list ${title} eta,xi,phi figure_phi_xi_eta &
#srun -n1 $sem_utils/utils/xsection/plot_slice_sphere_alpha_beta_phi_xi_eta.sh nc slice_sphere_220km.list ${title} eta,xi,phi figure_phi_xi_eta

wait