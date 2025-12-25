#!/bin/bash

n=$(ls events/*/DATA/CMTSOLUTION | wc -l)

cat<<EOF
# vtk DataFile Version 2.0
Source and Receiver VTK file
ASCII
DATASET UNSTRUCTURED_GRID
POINTS      $n float
EOF

ls events/*/DATA/CMTSOLUTION | xargs grep "^[xyz](m)" |\
  awk '{printf "%e ",$2/6371000.0; if(NR%3==0) printf "\n"}'
