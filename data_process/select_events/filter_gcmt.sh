#!/bin/bash

gcmt_catalog=GCMT/GCMT_1976-NOW.evt

#====== 1. 5.5<=M<=7.0 and after year 2000
#
min_mag=5.5
max_mag=7.2
min_t0="2000-01-01T00:00:00"

awk -F"|" -v min_mag=${min_mag} -v max_mag=${max_mag} -v min_t0=${min_t0} \
  '$1!~/^#/&&$2>=min_t0&&$6>=min_mag&&$6<=max_mag' ${gcmt_catalog} > events_sel_1.lst

#====== 2. lies within mesh domain
#
python sem_utils/meshfem3d/sem_test_inside_mesh.py Par_file events_sel_1.lst -o events_sel_1_with_xi_eta.csv

awk -F, 'NR>1 && $NF!~/nan/ && $NF>0{print $1}' events_sel_1_with_xi_eta.csv > evid_sel_2

#====== 3. no other events in gcmt catalog occur within 1 hour before or 0.5 hour after the origin time

cat /dev/null > evid_sel_3

tmp_list=$(mktemp -p .)
for evid in $(cat evid_sel_2)
do
  echo ======= $evid

  # epoch time
  t0=$(grep "^${evid}|" $gcmt_catalog | TZ=UTC awk -F"|" '{gsub("[-:T]"," ",$2); print mktime($2)}')

  # echo $t0

  TZ=UTC awk -F"|" -v evid=${evid} -v t0=${t0} '$1!~/^#/&&$1!=evid{gsub("[-:T]"," ",$2); t1=mktime($2); dt=t1-t0; if (dt>-3600 && dt<1800) print $1, dt }' $gcmt_catalog > $tmp_list

  cat $tmp_list
  if [ ! -s $tmp_list ];
  then
    echo $evid >> evid_sel_3
  fi

done

rm $tmp_list

grep -f evid_sel_3 $gcmt_catalog > events_sel_3.lst
awk -F"|" '{gsub(/ /,"",$5); print $3,$4,$1"_D"$5"_M"$6}' events_sel_3.lst | ./seis_utils/utils/mkKML.sh > evid_sel_3.kml


#====== 4. for any two events within a distance of 1 degree and origin time closer than 6 months, keep the largest one

python sem_utils/data_process/select_events/decluster.py events_sel_3.lst events_sel_4.lst --min_dist_in_km 100 --min_time_diff_in_days 180

awk -F"|" '{gsub(/ /,"",$5); print $3,$4,$1"_D"$5"_M"$6}' events_sel_4.lst | ./seis_utils/utils/mkKML.sh > evid_sel_4.kml
