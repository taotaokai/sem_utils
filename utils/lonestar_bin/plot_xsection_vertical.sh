#!/bin/bash

# plot profiles created from xsem_slice_...

#====== command line args
control_file=${1:?must provide control_file}
evid=${2:?must provide evid}
xsection_list=${3:?must provide xsection_list}
model_names=${4:-vsv,vsh,vpv,vph}

#====== set parameters
# source parameters in control_file
source ${control_file}
# get full path
xsection_list=$(readlink -f $xsection_list)
# etopo1
etopo1_grd=/Users/taok/research/GeoData/topo/ETOPO1_Bed_c_gmt4.grd

#====== plot each xsection
# create working directory
tmp_dir=$(mktemp -d -p.)
cd $tmp_dir
  
grep -v "^#" $xsection_list |\
while read lat0 lon0 azimuth theta0 theta1 ntheta r0 r1 nr fname
do
    echo "# xsection: $fname ($lat0 $lon0 $azimuth $theta0 $theta1 $ntheta $r0 $r1 $nr)"

    # input data files
    nc_file=$iter_dir/$evid/xsection/${fname}.nc
    ref_nc_file=$wkdir/1D_REF/xsection/${fname}.nc
    # output figure
    ps=$iter_dir/$evid/xsection/${fname}.ps

    gmt set PS_MEDIA A4
    gmt set PS_PAGE_ORIENTATION portrait

    #------ plot Map and surface traction of the cross_section
    R=90/150/20/55
    J=L120/38/30/45/10c
    #J=M120/38/8c

    # make basemap 
    gmt psbasemap -Yf22c -Xc -R$R -J$J -BWSNE -Bag -K > $ps
    gmt makecpt -Cglobe > cpt
    gmt grdimage $etopo1_grd -R -J -Ccpt -O -K >> $ps
    gmt pscoast -R -J -Dl -A250 -N1/thinnest -O -K >> $ps
    # plot cross-section with markers
    seq $theta0 5 $theta1 > theta.list
    $sem_utils/bin/xmap_gcircle $lat0 $lon0 $azimuth theta.list \
        > theta_lat_lon.list 
    awk '{print $3,$2}' theta_lat_lon.list |\
        gmt psxy -R -J -Wthin,black -O -K >> $ps
    # mark starting point
    awk 'NR==1{print $3,$2}' theta_lat_lon.list |\
        gmt psxy -R -J -Sc0.2c -Gred -O -K >> $ps
    awk 'NR>1{print $3,$2}' theta_lat_lon.list |\
        gmt psxy -R -J -Sc0.2c -Wthin -O -K >> $ps

    #------ plot cross-sections of vpv,vph,vsv,vsh
    R=$theta0/$theta1/$r0/$r1
    J=Pa12c/0z
    
    gmt grdmath -R$R -I1/1 Y = rr.nc
    echo 410 660 | awk '{printf "%f C\n%f C\n",6371-$1,6371-$2}' > c.txt
    
    Y=1
    for model_tag in ${model_names//,/ }
    do
        # calculate relative perturbation referred to 1D_REF
        gmt grdreformat $nc_file?$model_tag grd
        gmt grdreformat $ref_nc_file?$model_tag ref_grd
        gmt grdmath grd ref_grd SUB ref_grd DIV 100.0 MUL = dgrd
        # plot cross-section
        gmt grd2cpt dgrd -Cseis -R$R -Z -T= -D -L-5/5 > cpt
        gmt grdimage dgrd -Yf${Y}c -Xc -R$R -J$J -Ccpt \
            -BSWne -Bxa10f5 -Bya200f100 -O -K >> $ps
        # plot colorbar
        gmt psscale -D13c/2c/3.5c/0.2c -Ccpt \
            -B1:"dln($model_tag)":/:"(%)": -O -K >> $ps
        # mark starting point
        awk 'NR==1{printf "%f 6371\n",$1}' theta.list |\
            gmt psxy -R -J -Sc0.2c -Gred -N -O -K >> $ps
        awk 'NR>1{printf "%f 6371\n",$1}' theta.list |\
            gmt psxy -R -J -Sc0.2c -Wthin -N -O -K >> $ps
        # mark 410/660-km
        gmt grdcontour rr.nc -J$J -Cc.txt -O -K >> $ps
    
        # Y-shift
        Y=$((Y + 5))
    done

done # xsection_list

#====== clean
rm -rf $tmp_dir

#END