#!/bin/bash

# perl gmt2vtk.pl zhangpz_block.txt > zhangpz_block.vtk
# #perl gmt2vtk.pl zhangpz_pb.txt > zhangpz_pb.vtk
# perl gmt2vtk.pl pb.txt > pb.vtk
#

#R=100/160/15/60
R=-70/70/0/80

gmt pscoast -N1 -Dl -R$R -M > boundary.txt
perl gmt2vtk.pl boundary.txt > boundary.vtk

gmt pscoast -W1 -Dl -R$R -M > shoreline.txt
perl gmt2vtk.pl shoreline.txt > shoreline.vtk
