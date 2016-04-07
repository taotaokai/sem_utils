#!/bin/bash

# run xspecfem3D

wkdir=$(pwd)

data_dir=${1:?"[arg] need data_dir(for Par_file,CMTSOLUTION,STATION..)"}
sem_dir=${2:?"[arg] need sem_dir(for bin/xspecfem3D)"}
work_flow=${3:?"[arg] need work_flow(e.g. syn,adj,dstf,dxs,dmt,hess)"}

par_file=$data_dir/Par_file
adj_type=$(grep ^SIMULATION_TYPE $par_file | sed "s/[^=]*=[ ]*\([^ ]*\).*/\1/")
nxi=$(grep ^NPROC_XI $par_file | awk '{print $NF}')
neta=$(grep ^NPROC_ETA $par_file | awk '{print $NF}')
nproc=$(echo "$nxi $neta" | awk '{print $1*$2}')

function setup_dirs{
  out_dir=$1
  adj_type=$2
  # prepare output dir
  if [ -d "$out_dir" ]; then
    rm -rf $out_dir
  fi
  mkdir $out_dir
  rm OUTPUT_FILES
  ln -sf $out_dir OUTPUT_FILES
  cp DATABASES_MPI/addressing.txt OUTPUT_FILES
  # backup parameter files
  cp -L $data_dir/Par_file OUTPUT_FILES
  if [ x$adj_type == x1 ]; then
    cp -L $data_dir/STATIONS OUTPUT_FILES
  else
    cp -L $data_dir/STATIONS_ADJOINT OUTPUT_FILES
  fi
  cp -L $data_dir/CMTSOLUTION OUTPUT_FILES
}


for wk in $(sed "s/,/ /" $work_flow)
do

  if [ x"$wk" == xsyn ]
  then
    setup_dirs output_green 1
    ln -s $datadir
    echo mpirun -np $nproc $sem_dir/bin/xspecfem3D
    mpirun -np $nproc $sem_dir/bin/xspecfem3D
  fi

  if [ x"$wk" == xadj ]
  then
    setup_dirs output_srcfrechet
    echo mpirun -np $nproc $sem_dir/bin/xspecfem3D
    mpirun -np $nproc $sem_dir/bin/xspecfem3D
    # read in src_frechet file
python - <<EOF
misfit = Misfit()
misfit.load(filename="$misfit_dir/misfit.pkl")
misfit.read_srcfrechet(filename="output_srcfrechet/src_frechet.0001", \
  update=True)
misfit.save(filename='%s/misfit.pkl' % (misfit_dir))
EOF

  fi

  if [ x"$wk" == xdstf ]
  then
    echo mpirun -np $nproc $sem_dir/bin/xspecfem3D
    mpirun -np $nproc $sem_dir/bin/xspecfem3D
  fi


done

# run
echo mpirun -np $nproc $sem_dir/bin/xspecfem3D
mpirun -np $nproc $sem_dir/bin/xspecfem3D