#!/bin/bash

# make a list of spherical slices for xsem_slice_sphere

#---- 
sem_utils=./sem_utils

#---- mesh center
mesh_lat0=38.5
mesh_lon0=118

#---- xsem_slice_sphere:
lat0=15
lat1=60
nlat=451
lon0=90
lon1=150
nlon=601

#---- depths: 20 - 1000km 
seq 20 20 1000 > depth.list

#----
out_file="slice_sphere.list"

#---- make profile list
echo "# spherical slices:" > $out_file

cat depth.list |\
while read depth
do
    printf "%6.2f %6.2f %4d %7.2f %7.2f %4d %6.1f D_%06.1f\n" \
        $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth $depth
done >> $out_file
