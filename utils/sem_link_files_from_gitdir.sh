#!/bin/bash

gitdir=../specfem3d_globe
wkdir=$(pwd)


# clean 
find . -type l | xargs rm

# copy only necessary files 
ln -s $gitdir/config* ./
ln -s $gitdir/flags.guess ./
ln -s $gitdir/install-sh ./
ln -s $gitdir/Makefile.in ./
ln -s $gitdir/setup ./

cp -r $gitdir/src ./


# DATA
mkdir DATA

cd DATA
ln -s ../$gitdir/DATA/* ./
rm CMTSOLUTION STATIONS Par_file

mkdir control_files
cd control_files
cp ../../$gitdir/DATA/CMTSOLUTION ./
cp ../../$gitdir/DATA/STATIONS ./
cp ../../$gitdir/DATA/Par_file ./
