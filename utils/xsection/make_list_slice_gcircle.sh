#!/bin/bash

# get center points for x-sections in great circles

#---- 
sem_utils=${1:?[arg]need sem_utils directory for bin/xmap_gcircle}
out_file=${2:?[arg]need output file (e.g. slice_gcircle.list)}

#---- mesh center
mesh_lat0=38.5
mesh_lon0=118

#---- xsem_slice_gcircle:
# gcircle angular range, number of interpolations 

# radius range, number of interpolations
r0=5371
r1=6371
nr=101

flag_ellipticity=1

#----
#out_file="slice_gcircle.list"

#====== create lists of profile centers

##---- E-W gcircles:
#az=0
#seq -22 2 22 > theta.list
#$sem_utils/bin/xmap_gcircle \
#    $mesh_lat0 $mesh_lon0 \
#    $az theta.list > EW_center.list
#
##---- N-S gcircles:
#az=90
#seq -22 2 22 > theta.list
#$sem_utils/bin/xmap_gcircle \
#    $mesh_lat0 $mesh_lon0 \
#    $az theta.list > NS_center.list

#---- trench parallel gcircles:
az=116.5
seq -15 1 20 > theta.list
$sem_utils/bin/xmap_gcircle \
    $mesh_lat0 $mesh_lon0 \
    $az 1 theta.list > trench_parallel_center.list

#---- trench normal gcircles:
az=26.5
seq -22 1 21 > theta.list
$sem_utils/bin/xmap_gcircle \
    $mesh_lat0 $mesh_lon0 \
    $az 1 theta.list > trench_normal_center.list

#====== make profile list

## E-W gcircles
#echo "# E-W gcircles:" > $out_file
#
#az=90; theta0=-22; theta1=22; ntheta=441
#cat EW_center.list |\
#while read theta lat lon
#do
#    printf "%6.2f %7.2f %05.1f %6.2f %6.2f %4d %7.2f %7.2f %4d EW_%4.1f_%05.1f\n" \
#        $lat $lon $az $theta0 $theta1 $ntheta $r0 $r1 $nr $lat $lon
#done >> $out_file
#
## N-S gcircles
#echo "# N-S gcircles:" >> $out_file
#
#az=0; theta0=-22; theta1=22; ntheta=441
#cat NS_center.list |\
#while read theta lat lon
#do
#    printf "%6.2f %7.2f %05.1f %6.2f %6.2f %4d %7.2f %7.2f %4d NS_%4.1f_%05.1f\n" \
#        $lat $lon $az $theta0 $theta1 $ntheta $r0 $r1 $nr $lat $lon
#done >> $out_file

# trench parallel gcircles
echo "#TK trench parallel gcircles:" >> $out_file

az=26.5; theta0=-20; theta1=20; ntheta=401
cat trench_parallel_center.list |\
while read theta lat lon
do
    printf "%6.2f %7.2f %05.1f %6.2f %6.2f %4d %7.2f %7.2f %4d %d trench_parallel_lon%05.1f\n" \
        $lat $lon $az $theta0 $theta1 $ntheta $r0 $r1 $nr $flag_ellipticity $lon
done >> $out_file

# trench normal gcircles
echo "# trench normal gcircles:" >> $out_file

az=116.5; theta0=-15; theta1=20; ntheta=351
cat trench_normal_center.list |\
while read theta lat lon
do
    printf "%6.2f %7.2f %05.1f %6.2f %6.2f %4d %7.2f %7.2f %4d %d trench_normal_lat%04.1f\n" \
        $lat $lon $az $theta0 $theta1 $ntheta $r0 $r1 $nr $flag_ellipticity $lat
done >> $out_file
