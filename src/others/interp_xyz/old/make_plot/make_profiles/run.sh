#!/bin/bash

rm figure/*

profilename=section

matlab -nodesktop -nosplash <<EOF
plot_profile('$profilename','vs',200);
plot_profile('$profilename','vp',100);
exit
EOF

#plot_profile('$profilename','rho',200);
#plot_profile_relative('$profilename','vsh',200);
#plot_profile_relative('$profilename','vsv',200);
#plot_profile_relative('$profilename','vph',100);
#plot_profile_relative('$profilename','vpv',100);
#plot_profile_relative('$profilename','rho',200);
#plot_profile_relative('$profilename','eta',200);
#
