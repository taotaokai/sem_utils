#!/bin/bash

# print out weighted average cc of each event at the specified iteration steps 

wkdir=$(pwd)
event_list=${1:?[arg]need event_list}
iter_num=${2:?[arg]need iter_num e.g. 0,1,2,3}

niter=$(echo ${iter_num//,/ } | wc -w)

rm *.pdf *.list
rm iter??
ln -s ../iter?? .

#====== print wcc for all data
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
    END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no measurements
awk 'NF>1' list > all.list

./plot_wcc.py $niter all.list wcc_all.pdf

#====== print wcc for S body wave on transverse component
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2!~/surface/&&$2~/_T/&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
      END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no measurements
awk 'NF>1' list > body_wave_T.list

./plot_wcc.py $niter body_wave_T.list wcc_body_wave_T.pdf

#====== print wcc for P body wave on vertical component
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2!~/surface/&&($2~/P_Z/)&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
      END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no measurements
awk 'NF>1' list > body_wave_P-Z.list

./plot_wcc.py $niter body_wave_P-Z.list wcc_body_wave_P-Z.pdf

#====== print wcc for P body wave on radial component
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2!~/surface/&&($2~/P_R/)&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
      END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no measurements
awk 'NF>1' list > body_wave_P-R.list

./plot_wcc.py $niter body_wave_P-R.list wcc_body_wave_P-R.pdf

#====== print wcc for S body wave on vertical component
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2!~/surface/&&($2~/S_Z/)&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
      END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no measurements
awk 'NF>1' list > body_wave_S-Z.list

./plot_wcc.py $niter body_wave_S-Z.list wcc_body_wave_S-Z.pdf

#====== print wcc for S body wave on radial component
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2!~/surface/&&($2~/S_R/)&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
      END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no measurements
awk 'NF>1' list > body_wave_S-R.list

./plot_wcc.py $niter body_wave_S-R.list wcc_body_wave_S-R.pdf


#====== print wcc for surface_T
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2~/surface_T/&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
      END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no surface wave measurements
awk 'NF>1' list > surface_T.list

./plot_wcc.py $niter surface_T.list wcc_surface_T.pdf

#====== print wcc for surface_R
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2~/surface_R/&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
    END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no surface wave measurements
awk 'NF>1' list > surface_R.list

./plot_wcc.py $niter surface_R.list wcc_surface_R.pdf

#====== print wcc for surface_Z
awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list |\
while read event_id
do

  echo -n "$event_id"

  for iter in ${iter_num//,/ }
  do
    iter_dir=$(echo $iter | awk '{printf "iter%02d",$1}')
    misfit_file=$wkdir/$iter_dir/misfit_hist/${event_id}.txt

    # wcc
    awk '$1!~/#/&&$2~/surface_Z/&&$6<=10&&$6>=-10{n+=1; w+=$3; cc+=$3*$4; dt+=$3*$6} 
    END{if (n>0) printf " | %4d %6.1f %8.6f %+6.3f ", n, w, cc/w, dt/w}' $misfit_file
  done

  echo

done > list 

# delete events with no surface wave measurements
awk 'NF>1' list > surface_Z.list

./plot_wcc.py $niter surface_Z.list wcc_surface_Z.pdf

#====== merger figures
pdf_merge.sh wcc_hist.pdf wcc_*.pdf