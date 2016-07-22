#!/bin/bash

# make a list of spherical slices for xsem_slice_sphere

out_file=${1:?[arg]need output list (e.g. slice_sphere.list)}

#---- geographic range and mesh size
lat0=15
lat1=60
nlat=451
lon0=90
lon1=150
nlon=601

#---- depths(km): 5, 10, 20:20:1000 
seq 5 5 10 > depth.list
seq 20 20 900 >> depth.list

#---- make list
echo "# spherical slices:" > $out_file

cat depth.list |\
while read depth
do
    printf "%6.2f %6.2f %4d %7.2f %7.2f %4d %6.1f Depth_%03d\n" \
        $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth $depth
done >> $out_file
