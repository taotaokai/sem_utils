#!/bin/bash

# run xspecfem3D
datefmt='-u +%Y-%m-%dT%H:%M:%S'
wkdir=$(pwd)

sem_dir=${1:?"[arg] need sem_dir(for bin/xspecfem3D)"}
work_flow=${2:?"[arg] need work_flow(e.g. syn,adj,dstf,dxs,dmt,hess)"}

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

  if [ x"$wk" == xsyn ]; then
    echo ====== $wk
    out_dir=output_green
    rm -rf $out_dir OUTPUT_FILES
    mkdir $out_dir
    ln -sf $out_dir OUTPUT_FILES
    cp DATABASES_MPI/addressing.txt OUTPUT_FILES
    cp -L DATA/Par_file OUTPUT_FILES
    cp -L DATA/CMTSOLUTION OUTPUT_FILES
    cp -L DATA/STATIONS OUTPUT_FILES
    sed -i "/^[\s]*SIMULATION_TYPE/s/=.*/= 1/" $par_file
    echo [$(date $datefmt)] mpirun -np $nproc $sem_dir/bin/xspecfem3D
    mpirun -np $nproc $sem_dir/bin/xspecfem3D
    echo [$(date $datefmt)] done.
    sleep 1m
  fi

#  if [ x"$wk" == xadj ]
#  then
#    setup_dirs output_srcfrechet
#    echo mpirun -np $nproc $sem_dir/bin/xspecfem3D
#    mpirun -np $nproc $sem_dir/bin/xspecfem3D
#    # read in src_frechet file
#python - <<EOF
#misfit = Misfit()
#misfit.load(filename="$misfit_dir/misfit.pkl")
#misfit.read_srcfrechet(filename="output_srcfrechet/src_frechet.0001", \
#  update=True)
#misfit.save(filename='%s/misfit.pkl' % (misfit_dir))
#EOF
#
#  fi
#
#  if [ x"$wk" == xdstf ]
#  then
#    echo mpirun -np $nproc $sem_dir/bin/xspecfem3D
#    mpirun -np $nproc $sem_dir/bin/xspecfem3D
#  fi


done