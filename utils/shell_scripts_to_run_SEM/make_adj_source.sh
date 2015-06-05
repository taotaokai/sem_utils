#!/bin/bash

#calculate the adjoint sources
#echo "# $(date)"

# read Par_file
source Par_file 

# check if mehser finished successfully
#mesher_status=$(awk -F"=" '/done/{print $2}' $iter_dir/mesher.status)
#if [ -z "$mesher_status" ] || [ "$mesher_status" -ne 1 ]
#then
#	echo "ERROR: Run run_mesher.sh first!"
#	exit
#fi

#
for evnm in $(cat $iter_dir/DATA/EVENTS)
do
    evnm_dir=$iter_dir/EVENTS/$evnm
    echo $evnm_dir

    # check if forward finished successfully
    #forward_status=$(awk -F"=" '/done/{print $2}' $evnm_dir/forward.status)
    #if [ -z "$forward_status" ] || [ "$forward_status" -ne 1 ]
    #then
   	#    echo ERROR: Run forward simutlation first!
    #    exit
    #fi
    
    # set sac header
    echo "set sac header"
    cd $evnm_dir/OUTPUT_forward
sac<<EOF
r *.sac
wh
q
EOF

    # make adjoint sources
    echo "measure adjoint sources"
    cd $evnm_dir
    $USER_BIN/measure_adj_src.py DATA/STATIONS \
                                 $event_dir/$evnm/DATA/PHASES_ADJOINT \
                                 $event_dir/$evnm/DISPLACEMENT \
                                 OUTPUT_forward \
                                 ADJOINT_SOURCES \
                                 $freqmin $freqmax > ADJOINT_SOURCES/adj.out
    # sort adj.out with cc
    sort -k14 -n ADJOINT_SOURCES/adj.out > tmp
    mv tmp ADJOINT_SOURCES/adj.out

    # generate STATION_ADJOINT
    #choose only stations of cc >= cutoff(this value is now chosen arbitarily)
    awk '$14>c{printf "%5s %5s %8s %9s %8s   0\n",$2,$1,$3,$4,$5}' c=$cc_cutoff \
        ADJOINT_SOURCES/adj.out >  DATA/STATIONS_ADJOINT
    cp DATA/STATIONS_ADJOINT ADJOINT_SOURCES/
done

echo "DONE"

#END
