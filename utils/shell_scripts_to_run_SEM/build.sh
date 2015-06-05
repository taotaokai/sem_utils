#!/bin/bash

# configure and compile all the SEM programs

# read Par_file
source Par_file

# sets up directory structure in current example directoy
echo "configure and build"
echo $(date)
echo "SEM_dir=$SEM_dir"
echo "base_dir=$base_dir"

# configure
cp $config_dir/setup/*.h.in $SEM_dir/setup
cd $SEM_dir
echo "configure ..."
./configure FC=gfortran MPIFC=mpif90 2>&1 | tee $config_dir/configure.log

# build src/
cp $config_dir/DATA/Par_file $SEM_dir/DATA/
echo "make ..."
make clean 2>&1 | tee $config_dir/make.log
make all 2>&1 | tee -a $config_dir/make.log
make tomography 2>&1 | tee -a $config_dir/make.log

# utils
if [ -d $SEM_dir/utils/Visualization/VTK_ParaView/mesh2vtu ]
then
  cd $SEM_dir/utils/Visualization/VTK_ParaView/mesh2vtu
  make clean
  make
fi

# user contributed
if [ -d $SEM_dir/src/user_contrib/mask_source ]
then
  cd $SEM_dir/src/user_contrib/mask_source
  make clean
  make
fi

# backup parameter files
mkdir $config_dir/OUTPUT_FILES
cd $SEM_dir
cp setup/*.h $config_dir/OUTPUT_FILES
cp OUTPUT_FILES/values_from_mesher.h $config_dir/OUTPUT_FILES

# copy executables
mkdir $base_dir/SEM_bin
rm -rf $base_dir/SEM_bin/*
cp $SEM_dir/bin/x* $base_dir/SEM_bin/
cp $SEM_dir/utils/Visualization/VTK_ParaView/mesh2vtu/mesh2vtu $base_dir/SEM_bin/
cp $SEM_dir/src/user_contrib/mask_source/xmask_source $base_dir/SEM_bin/

# link DATA files
cd $config_dir/DATA/
#find . -xtype l | awk '{printf "rm %s\n",$1}' | bash  # remove broken links
find . -type l | xargs rm # remove all symlinks
ln -s $SEM_dir/DATA/crust2.0
ln -s $SEM_dir/DATA/s362ani
ln -s $SEM_dir/DATA/QRFSI12
ln -s $SEM_dir/DATA/topo_bathy

echo $(date)
echo "DONE."
#END
