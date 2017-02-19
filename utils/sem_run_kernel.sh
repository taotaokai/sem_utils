#!/bin/bash

# setup and generate the batch script to run xspecfem3D for adjoint simulation (kernel)

#proc_per_node=24
proc_per_node=20

#====== command line args
run_dir=${1:?[arg]need run_dir(for all output)}
syn_dir=${2:?[arg]need syn_dir(for DATA/Par_file,STATIONS,CMTSOLUTION; DATABASES_MPI/*, save_forward=.true.)}
adj_dir=${3:?[arg]need adj_dir(for adjoint source files *.adj)}
sem_dir=${4:?[arg]need sem_dir(for code, DATA/*)}
mpiexec=${5:?[arg]need mpiexec (e.g. ibrun or mpirun -np 144)}

if [ -d "$run_dir" ]
then
    echo "[WARN] run_dir($run_dir) exists, delete!"
    rm -rf $run_dir
fi
mkdir $run_dir

if [ ! -d "$syn_dir" ]
then
    echo "[ERROR] syn_dir($syn_dir) does NOT exist!"
    exit 1
elif [ ! -f "$syn_dir/DATA/Par_file" ]
then
    echo "[ERROR] $syn_dir/DATA/Par_file does NOT exist!"
    exit 1
elif [ ! -f "$syn_dir/DATA/CMTSOLUTION" ]
then
    echo "[ERROR] $syn_dir/DATA/CMTSOLUTION does NOT exist!"
    exit 1
elif [ ! -f "$syn_dir/DATA/STATIONS" ]
then
    echo "[ERROR] $syn_dir/DATA/STATIONS does NOT exist!"
    exit 1
fi

if [ ! -d "$adj_dir" ]
then
    echo "[ERROR] adj_dir($adj_dir) does NOT exit!"
    exit 1
fi

run_dir=$(readlink -f $run_dir)
syn_dir=$(readlink -f $syn_dir)
adj_dir=$(readlink -f $adj_dir)
sem_dir=$(readlink -f $sem_dir)

#====== setup run_dir
cd $run_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

# link data files
cd $run_dir/DATA

rm -rf Par_file STATIONS CMTSOLUTION
cp -L $syn_dir/DATA/Par_file .
cp -L $syn_dir/DATA/CMTSOLUTION .
cp -L $syn_dir/DATA/STATIONS .

sed -i "/^SIMULATION_TYPE/s/=.*/= 3/" Par_file
sed -i "/^SAVE_FORWARD/s/=.*/= .false./" Par_file
sed -i "/^ANISOTROPIC_KL/s/=.*/= .false./" Par_file
sed -i "/^SAVE_TRANSVERSE_KL_ONLY/s/=.*/= .false./" Par_file
sed -i "/^APPROXIMATE_HESS_KL/s/=.*/= .false./" Par_file
#sed -i "/^[\s]*SAVE_MESH_FILES/s/=.*/= .false./" Par_file
#sed -i "/^[\s]*MODEL/s/=.*/= GLL/" Par_file

# make STATIONS_ADJOINT
cd $adj_dir
#ls *Z.adj | sed 's/..Z\.adj$//' |\
#  awk -F"." '{printf "%s[ ]*%s.%s[ ]\n",$1,$2,$3}' > $run_dir/DATA/grep_pattern
ls *Z.adj | sed 's/\...Z\.adj$//' |\
  awk -F"." '{st=index($0,"."); printf "%s[ ]*%s\n",$1,substr($0,st+1)}' > $run_dir/DATA/grep_pattern
grep -f $run_dir/DATA/grep_pattern $syn_dir/DATA/STATIONS \
  > $run_dir/DATA/STATIONS_ADJOINT

# backup DATA/* into OUTPUT_FILES/
cd $run_dir/DATA
cp -L Par_file CMTSOLUTION STATIONS STATIONS_ADJOINT $run_dir/OUTPUT_FILES/

# link adjoint source dir
ln -s $adj_dir $run_dir/SEM

# link mesh database
cd $run_dir/DATABASES_MPI
ln -s $syn_dir/DATABASES_MPI/*.bin .

# OUTPUT_FILES
cp $syn_dir/OUTPUT_FILES/addressing.txt $run_dir/OUTPUT_FILES
cp $syn_dir/OUTPUT_FILES/addressing.txt $run_dir/DATABASES_MPI

# generate sbatch job file
nproc_xi=$(grep NPROC_XI $run_dir/DATA/Par_file | awk '{print $NF}')
nproc_eta=$(grep NPROC_ETA $run_dir/DATA/Par_file | awk '{print $NF}')
nproc=$(echo "$nproc_xi * $nproc_eta" | bc -l)
nnode=$(echo "$nproc $proc_per_node" | awk '{a=$1/$2}END{print (a==int(a))?a:int(a)+1}')
 
cat <<EOF > $run_dir/kernel.job
#!/bin/bash
#SBATCH -J kernel
#SBATCH -o $run_dir/kernel.job.o%j
#SBATCH -N $nnode
#SBATCH -n $nproc
#SBATCH -p normal
#SBATCH -t 01:00:00
#SBATCH --mail-user=kai.tao@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd $run_dir
${mpiexec} $sem_dir/bin/xspecfem3D

EOF
#END
