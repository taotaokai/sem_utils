#!/bin/bash

wkdir=$(pwd)

event_list=${1:?[arg]need event_list}

for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
do

  echo "====== $event_id"

  utils/ecef2lla.py CMTSOLUTION_initial/$event_id.cmt
  utils/ecef2lla.py CMTSOLUTION_updated/$event_id.cmt

  x0=$(grep 'x(m)' CMTSOLUTION_initial/$event_id.cmt | awk '{print $2}')
  y0=$(grep 'y(m)' CMTSOLUTION_initial/$event_id.cmt | awk '{print $2}')
  z0=$(grep 'z(m)' CMTSOLUTION_initial/$event_id.cmt | awk '{print $2}')

  x1=$(grep 'x(m)' CMTSOLUTION_updated/$event_id.cmt | awk '{print $2}')
  y1=$(grep 'y(m)' CMTSOLUTION_updated/$event_id.cmt | awk '{print $2}')
  z1=$(grep 'z(m)' CMTSOLUTION_updated/$event_id.cmt | awk '{print $2}')

  #echo "(($x1 - $x0)^2 + ($x1 - $x0)^2 + ($x1 - $x0)^2)^0.5/1000.0" | bc -l
  # displacement in source location (km)
  echo $x0 $x1 $y0 $y1 $z0 $z1 | awk '{print (($1-$2)^2+($3-$4)^2+($5-$6)^2)^0.5/1000.0}'

done > check_lla.out

#awk 'BEGIN{RS="======"}{printf "%s %8.5f %8.5f %8.5f\n", $1, $5-$2, $6-$3, $7-$4}' check_lla.out > diff_lla.out

# dlat, dlon(deg), ddep(km)
#awk 'BEGIN{RS="======"; pi=atan2(0,-1); deg2rad=pi/180} {printf "%s %8.5f %8.5f %8.5f\n", $1, $5-$2, $6-$3, $7-$4, $8}' check_lla.out > diff_lla.out
echo "# event_id dlat(deg) dlon(deg) ddep(km) dloc(km)" > diff_lla.out
awk 'BEGIN{RS="======";} {printf "%s %8.5f %8.5f %8.5f %8.5f\n", $1, $5-$2, $6-$3, $7-$4, $8}' check_lla.out | sed "1d" | sort -n -k5 >> diff_lla.out