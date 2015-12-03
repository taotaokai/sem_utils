#!/bin/bash

# plot profiles created from xsem_slice_...

# A4 paper: H29.7cm x W21cm

#====== command line args
control_file=${1:?[arg] need control_file}
event_id=${2:?[arg] need event_id}
slice_list=${3:?[arg] need slice_list}
model_names=${4:-vsv,vsh,vpv,vph}
ref_dir=${5:-1D_REF} # reference model

# load parameters in control_file
source ${control_file}

# get full path
slice_list=$(readlink -f $slice_list)

# etopo1
etopo1_grd=$base_dir/topo/ETOPO1_Bed_c_gmt4.grd

#====== plot each xsection

# create working directory
export GMT_TMPDIR=$(mktemp -d -p"$(pwd)")
echo $GMT_TMPDIR
cd $GMT_TMPDIR

# get number of xsections
xsection_num=$(echo ${model_names//,/ } | awk '{print NF}')
# height of each xsection (A4 paper: 29.7cm x 21cm)
xsection_yshift=$(echo "scale=4; 28/${xsection_num}" | bc -l)
xsection_height=$(echo "scale=4; ${xsection_yshift} * 0.8" | bc -l)
cbar_length=$(echo "scale=4; ${xsection_height} * 0.7" | bc -l)
cbar_ypos=$(echo "scale=4; ${xsection_height} * 0.5" | bc -l)

grep -v "^#" $slice_list |\
while read lat0 lat1 nlat lon0 lon1 nlon depth fname
do
    echo
    echo "#====== $fname: $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth"
    echo

    # input data files
    nc_file=$iter_dir/$event_id/xsection/${fname}.nc
    ref_nc_file=$base_dir/${ref_dir}/xsection/${fname}.nc
    # output figure
    ps=$iter_dir/$event_id/xsection/${fname}_STW105.ps

    gmt set \
        FONT_ANNOT_PRIMARY +12p,Times-Roman \
        FONT_TITLE 14p,Times-Bold \
        PS_PAGE_ORIENTATION portrait \
        MAP_FRAME_TYPE plain \
        PS_MEDIA A4

    #------ plot spherical slices for each depth

    # title
    echo "Depth $depth (km)" | \
        gmt pstext -Xc0c -Yc14c -R-5/5/-5/5 -JX10 \
        -F+cCM+f17,Times-Bold,+jCM -P -K > $ps

    R=$lon0/$lon1/$lat0/$lat1
    proj_lon0=$(echo "scale=4; ($lon0 + $lon1)/2.0" | bc -l)
    proj_lat0=$(echo "scale=4; ($lat0 + $lat1)/2.0" | bc -l)
    proj_lat1=$(echo "scale=4; $lat0 + ($lat1 - $lat0)/4.0" | bc -l)
    proj_lat2=$(echo "scale=4; $lat1 - ($lat1 - $lat0)/4.0" | bc -l)
    J=L${proj_lon0}/${proj_lat0}/${proj_lat1}/${proj_lat2}/${xsection_height}ch

    # plot each model parameter
    Y=1
    for model_tag in ${model_names//,/ }
    do
        echo "# $model_tag"
        # calculate relative perturbation referred to 1D_REF
        gmt grdreformat $nc_file?$model_tag grd
        gmt grdreformat $ref_nc_file?$model_tag ref_grd
        gmt grdmath grd ref_grd SUB ref_grd DIV 100.0 MUL = dgrd
        # plot cross-section
        gmt makecpt -Cjet -D -T-6/6/0.1 -Z -I > cpt
        #gmt makecpt -Cseis -D -T-5/5/0.1 -Z > cpt
        #gmt grd2cpt dgrd -Cjet -R$R -Z -D -I -L-6/6> cpt
        #gmt grd2cpt dgrd -Cseis -R$R -Z -T= -D > cpt
        gmt psbasemap -Yf${Y}c -Xc -R$R -J$J \
            -BWSne -Bxaf -Byaf -O -K >> $ps
        gmt grdimage dgrd -R -J -Ccpt -nb -O -K >> $ps
        gmt pscoast -R -J -Dl -A250 -N1/thinnest -W1/thinnest -O -K >> $ps
        # plot colorbar
        gmt psscale -D0c/${cbar_ypos}c/${cbar_length}c/0.15c -Ccpt \
            -Xf16c -Yf${Y}c \
            -Bxaf -By+l"d${model_tag}(%)" -E -O -K >> $ps

        # Y-shift
        Y=$(echo "scale=4; $Y + ${xsection_yshift}" | bc -l)
    done

    #------ covert .ps to .pdf file
    echo "#-- convert ps to pdf"
    ps2pdf $ps $iter_dir/$event_id/xsection/${fname}_STW105.pdf

done # xsection_list

#====== clean
rm -rf $GMT_TMPDIR

#END