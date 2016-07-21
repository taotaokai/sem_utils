#!/bin/bash

# plot profiles created from xsem_slice_...

#====== command line args
control_file=${1:?[arg] need control_file}
event_id=${2:?[arg] need event_id}
slice_list=${3:?[arg] need slice_list}
model_names=${4:-vsv,vsh,vpv,vph}

# load parameters in control_file
source ${control_file}

# get full path
slice_list=$(readlink -f $slice_list)

# etopo1
etopo1_grd=/Users/taok/research/GeoData/topo/ETOPO1_Bed_c_gmt4.grd

#====== plot each xsection

# create working directory
export GMT_TMPDIR=$(mktemp -d -p"$(pwd)")
echo $GMT_TMPDIR
cd $GMT_TMPDIR

# get number of xsections
xsection_num=$(echo ${model_names//,/ } | awk '{print NF}')
# height of each xsection (A4 paper: 29.7cm x 21cm)
xsection_yshift=$(echo "scale=4; 21/${xsection_num}" | bc -l)
xsection_height=$(echo "scale=4; ${xsection_yshift} * 0.8" | bc -l)
cbar_length=$(echo "scale=4; ${xsection_height} * 0.7" | bc -l)
cbar_ypos=$(echo "scale=4; ${xsection_height} * 0.45" | bc -l)

grep -v "^#" $slice_list |\
while read lat0 lon0 azimuth theta0 theta1 ntheta r0 r1 nr fname
do
    echo
    echo "#====== $fname: $lat0 $lon0 $azimuth $theta0 $theta1 $ntheta $r0 $r1 $nr"
    echo

    # input data files
    nc_file=$iter_dir/$event_id/xsection/${fname}.nc
    # output figure
    ps=$iter_dir/$event_id/xsection/${fname}.ps

    gmt set \
        FONT_ANNOT_PRIMARY +14p,Times-Roman \
        FONT_TITLE 14p,Times-Bold \
        PS_PAGE_ORIENTATION portrait \
        PS_MEDIA A4

    #------ plot Map and surface traction of the cross_section
    echo "#-- plot Map"

    R=90/150/20/55
    #J=L120/38/30/45/10c+
    J=T120/40/6ch
    #J=M120/38/8c
    
    # xsection lat,lon
    gmt convert ${nc_file}?theta > theta.list
    $sem_utils/bin/xmap_gcircle $lat0 $lon0 ${azimuth} theta.list > gcircle.list

    # xsection marker on map
    seq $theta0 5 $theta1 > theta.list
    $sem_utils/bin/xmap_gcircle $lat0 $lon0 ${azimuth} theta.list \
        > theta_lat_lon.list
    cat theta_lat_lon.list |\
    awk 'NR==1{printf "> -Sc0.2c -Gred -Wthin\n %s %s\n", $3, $2}; \
         NR==2{printf "> -Sc0.2c -Gwhite -Wthin\n %s %s\n", $3, $2};  \
         NR>2 {printf "%s %s\n", $3, $2}' > xsection_marker.xy

    # make basemap 
    gmt psbasemap -Yf22.5c -Xc -R$R -J$J -BWSNE -Bag -K > $ps
    gmt makecpt -Cglobe > cpt
    gmt grdimage $etopo1_grd -R -J -Ccpt -O -K >> $ps
    gmt pscoast -R -J -Dl -A250 -N1/thinnest -O -K >> $ps
    # gcircle
    awk '{print $3, $2}' gcircle.list | gmt psxy -R -J -Wthin -O -K >> $ps
    # xsection marker
    gmt psxy xsection_marker.xy -Sp -N -R -J -O -K >> $ps

    #------ plot cross-sections of vpv,vph,vsv,vsh
    echo "#-- plot xsection"

    # rotation angle for plotting slice
    angle=$(echo $theta0 $theta1 | awk '{print ($2-$1)/2.0 + $1}')

    R=$theta0/$theta1/$r0/$r1
    J=Pa${xsection_height}c/${angle}zh
    #J=Pa15c/${angle}z

    # lines for 410/660-km
    gmt grdmath -R$R -I1/1 Y = rr.nc
    echo 410 660 | awk '{printf "%f C\n%f C\n",6371-$1,6371-$2}' > c.txt

    # xsection marker
    cat theta_lat_lon.list |\
    awk 'NR==1{printf "> -Sc0.2c -Gred -Wthin\n %s 6371\n", $1}; \
         NR==2{printf "> -Sc0.2c -Gwhite -Wthin\n %s 6371\n", $1};  \
         NR>2 {printf "%s  6371\n", $1}' > xsection_marker.xy  

    # plot each model xsection
    Y=1
    for model_tag in ${model_names//,/ }
    do
        echo "# $model_tag"
        # plot cross-section
        gmt grdreformat ${nc_file}?${model_tag} grd
        #gmt makecpt -Cjet -D -T-6/6/0.1 -Z -I > cpt
        #gmt makecpt -Cseis -D -T-5/5/0.1 -Z > cpt
        gmt grd2cpt grd -Cjet -R$R -Z -D -I > cpt
        #gmt grd2cpt dgrd -Cseis -R$R -Z -T= -D > cpt
        gmt grdimage grd -Yf${Y}c -Xc -R$R -J$J -Ccpt -nb \
            -BWSne -Bxa10f5 -Bya200f100 -O -K >> $ps
        # xsection marker 
        gmt psxy xsection_marker.xy -Sp -N -R -J -O -K >> $ps
        # plot 410/660-km
        gmt grdcontour rr.nc -J$J -Cc.txt -O -K >> $ps
        # plot colorbar
        gmt psscale -D15c/${cbar_ypos}c/${cbar_length}c/0.15c -Ccpt \
            -Bxaf -By+l"${model_tag}" -E -O -K >> $ps

        # Y-shift
        Y=$(echo "scale=4; $Y + ${xsection_yshift}" | bc -l)

    done

    #------ covert .ps to .pdf file
    echo "#-- convert ps to pdf"
    ps2pdf $ps $iter_dir/$event_id/xsection/${fname}.pdf

done # slice_list

#====== clean
rm -rf $GMT_TMPDIR

#END