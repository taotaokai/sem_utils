#!/bin/bash

# run xspecfem3D
datefmt='-u +%Y-%m-%dT%H:%M:%S'
wkdir=$(pwd)
utils=$wkdir/utils

sem_dir=${1:?"[arg] need sem_dir(for bin/xspecfem3D)"}
work_flow=${2:?"[arg] need work_flow(e.g. syn,adj,der,hess)"}

if [ ! -d DATA ]; then
  echo "DATA/ must exist!"
  exit -1
fi
if [ ! -d DATABASES_MPI ]; then
  echo "DATABASES_MPI/ must exist!"
  exit -1
fi

par_file=DATA/Par_file
nxi=$(grep ^NPROC_XI $par_file | awk '{print $NF}')
neta=$(grep ^NPROC_ETA $par_file | awk '{print $NF}')
nproc=$(echo "$nxi $neta" | awk '{print $1*$2}')

for wk in ${work_flow//,/ }
do

  if [ x"$wk" == xsyn ]
  then
    echo ====== $wk
    out_dir=output_green
    rm -rf $out_dir OUTPUT_FILES
    mkdir $out_dir
    ln -sf $out_dir OUTPUT_FILES
    cd $wkdir/DATA
    cp CMTSOLUTION.init CMTSOLUTION.green
    sed -i "s/^tau(s):.*/tau(s):             0.0000000E+00/" CMTSOLUTION.green
    ln -sf CMTSOLUTION.green CMTSOLUTION
    cd $wkdir
    cp DATABASES_MPI/addressing.txt OUTPUT_FILES
    cp -L DATA/Par_file OUTPUT_FILES
    cp -L DATA/CMTSOLUTION OUTPUT_FILES
    cp -L DATA/STATIONS OUTPUT_FILES
    sed -i "/^[\s]*SIMULATION_TYPE/s/=.*/= 1/" $par_file
    echo [$(date $datefmt)] mpirun -np $nproc $sem_dir/bin/xspecfem3D
    mpirun -np $nproc $sem_dir/bin/xspecfem3D
    echo [$(date $datefmt)] done.
  fi

  if [ x"$wk" == xadj ]
  then
    echo ====== $wk
    # measure/plot misfit
    mkdir adj misfit
    ln -s adj SEM
    $utils/sac_mod.sh output_green "*.sac" "ch lcalda false; wh"
    $utils/measure_misfit.py
    $utils/plot_misfit.py
    # sem
    out_dir=output_srcfrechet
    rm -rf $out_dir OUTPUT_FILES
    mkdir $out_dir
    ln -sf $out_dir OUTPUT_FILES
    cd $wkdir/DATA
    rm CMTSOLUTION CMTSOLUTION.green
    cp CMTSOLUTION.init CMTSOLUTION.green
    sed -i "s/^tau(s):.*/tau(s):             0.0000000E+00/" CMTSOLUTION.green
    ln -sf CMTSOLUTION.green CMTSOLUTION
    ln -sf STATIONS STATIONS_ADJOINT
    cd $wkdir
    cp DATABASES_MPI/addressing.txt OUTPUT_FILES
    cp -L DATA/Par_file OUTPUT_FILES
    cp -L DATA/CMTSOLUTION OUTPUT_FILES
    cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
    sed -i "/^[\s]*SIMULATION_TYPE/s/=.*/= 2/" $par_file
    echo [$(date $datefmt)] mpirun -np $nproc $sem_dir/bin/xspecfem3D
    mpirun -np $nproc $sem_dir/bin/xspecfem3D
    echo [$(date $datefmt)] done.
  fi

  if [ x"$wk" == xder ]
  then
    echo "====== $wk (waveform der)"
    # make cmt of different source parameters
    $utils/make_cmt_der.py
    for der in dxs dmt
    do
      out_dir=output_$der
      rm -rf $out_dir OUTPUT_FILES
      mkdir $out_dir
      ln -sf $out_dir OUTPUT_FILES
      cp DATABASES_MPI/addressing.txt OUTPUT_FILES
      cd $wkdir/DATA
      ln -sf CMTSOLUTION.$der CMTSOLUTION
      cd $wkdir
      cp -L DATA/Par_file OUTPUT_FILES
      cp -L DATA/CMTSOLUTION OUTPUT_FILES
      cp -L DATA/STATIONS_ADJOINT OUTPUT_FILES
      sed -i "/^[\s]*SIMULATION_TYPE/s/=.*/= 1/" $par_file
      echo [$(date $datefmt)] mpirun -np $nproc $sem_dir/bin/xspecfem3D
      mpirun -np $nproc $sem_dir/bin/xspecfem3D
      echo [$(date $datefmt)] done.
    done
    # read in differential seismograms
    $utils/sac_mod.sh output_dxs "*.sac" "ch lcalda false; wh"
    $utils/sac_mod.sh output_dmt "*.sac" "ch lcalda false; wh"
    $utils/waveform_der.py
  fi

done