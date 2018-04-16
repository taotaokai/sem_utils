#!/bin/bash

perl gmt2vtk.pl zhangpz_block.txt > zhangpz_block.vtk 
#perl gmt2vtk.pl zhangpz_pb.txt > zhangpz_pb.vtk 
perl gmt2vtk.pl pb.txt > pb.vtk 

gmt pscoast -N1 -Dl -R100/160/15/60 -M > boundary.txt
perl gmt2vtk.pl country_border.txt > boundary.vtk

gmt pscoast -W1 -Dl -R100/160/15/60 -M > shoreline.txt
perl gmt2vtk.pl shoreline.txt > shoreline.vtk