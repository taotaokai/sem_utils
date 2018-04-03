#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np 

import matplotlib as mpl
#mpl.rcParams['font.family'] = "times new roman"
#mpl.rcParams['font.family'] = "times"
mpl.use("pdf")
import matplotlib.pyplot as plt

#====== input
input_list = str(sys.argv[1])
out_fig = str(sys.argv[2])

#====== read file list
with open(input_list, 'r') as f:
  file_list = [ l.split()[0] for l in f.readlines() if '#' not in l ]

#====== read in ref model
ref_file = "stage09.iter16.polyfit_um.linfit_mtz.vp_vs_voigt"
with open(ref_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_ref = np.array([float(l[1]) for l in lines])
vp_ref = np.array([float(l[2]) for l in lines])/1000.0
vs_ref = np.array([float(l[3]) for l in lines])/1000.0

FWEA18_1d_mean_file = "FWEA18_1d_mean.txt"
with open(FWEA18_1d_mean_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_FWEA18_1d_mean = np.array([float(l[0]) for l in lines])
vp_FWEA18_1d_mean = np.array([float(l[1]) for l in lines])
vs_FWEA18_1d_mean = np.array([float(l[2]) for l in lines])

stw105_file = "STW105.txt"
with open(stw105_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
r_stw105   = np.array([float(l[0]) for l in lines])/1000.
dep_stw105 = 6371. - r_stw105
vpv_stw105 = np.array([float(l[2]) for l in lines])/1000.0
vsv_stw105 = np.array([float(l[3]) for l in lines])/1000.0
vsh_stw105 = np.array([float(l[7]) for l in lines])/1000.0
vs_stw105 = ((2*vsv_stw105**2 + vsh_stw105**2)/3)**0.5

ak135_file = "AK135F.txt"
with open(ak135_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_ak135 = np.array([float(l[0]) for l in lines])
vp_ak135 = np.array([float(l[2]) for l in lines])
vs_ak135 = np.array([float(l[3]) for l in lines])

prem_file = "PREM_1s.txt"
with open(prem_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_prem = np.array([float(l[1]) for l in lines])
vpv_prem = np.array([float(l[3]) for l in lines])
vsv_prem = np.array([float(l[5]) for l in lines])
vsh_prem = np.array([float(l[6]) for l in lines])

vs_prem = ((2*vsv_prem**2 + vsh_prem**2)/3)**0.5

#====== read data files
with open(file_list[0], 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]

nmodel = len(file_list)
npts = len(lines)

r = np.zeros((npts, nmodel))
vsv = np.zeros((npts, nmodel))
vsh = np.zeros((npts, nmodel))
vpv = np.zeros((npts, nmodel))
vph = np.zeros((npts, nmodel))

imodel = 0
for data_file in file_list:
  with open(data_file, 'r') as f:
    lines = [ l.split() for l in f.readlines() if '#' not in l ]
  r[:, imodel]   = np.array([float(l[0]) for l in lines])
  vsv[:, imodel] = np.array([float(l[1]) for l in lines])
  vsh[:, imodel] = np.array([float(l[2]) for l in lines])
  vpv[:, imodel] = np.array([float(l[3]) for l in lines])
  vph[:, imodel] = np.array([float(l[4]) for l in lines])
  imodel += 1

vs = ((2*vsv**2 + vsh**2)/3)**0.5

#====== average model
vs_median = np.median(vs, axis=1)
#vp_median = np.median(vpv, axis=1)

vs_mean = np.mean(vsv, axis=1)
#vp_mean = np.mean(vpv, axis=1)

#====== plot Vsv
plt.figure()
#for imodel in range(nmodel):
#  plt.plot(vs[:, imodel], (1.0-r[:, imodel])*6371., 'k-', linewidth=0.2, alpha=0.1)

#line_median, = plt.plot(vs_median, (1.0-r[:,0])*6371., 'r-', linewidth=1.0)
line_mean, = plt.plot(vs_mean, (1.0-r[:,0])*6371., 'r-', linewidth=1.0)

with open('1d_mean_model.txt', 'w') as f:
  f.write("#1D mean model\n")
  f.write("#depth[km] vp[km/s] vs[km/s]\n")
  for i in range(npts):
    f.write("%12.6f  %10.5f  %10.5f\n"%((1.0-r[i,0])*6371.0, vp_mean[i], vs_mean[i]))

#plt.plot(vsv_mean, (1.0-r[:,0])*6371., 'b', linewidth=0.5)
#line_mean, = plt.plot(vsv_mean, (1.0-r[:,0])*6371., c='purple', ls='-', marker='.', lw=1.0)
#line_median, = plt.plot(vsv_median, (1.0-r[:,0])*6371., 'r-', linewidth=1.0)
#line_median, = plt.plot(kappav_median, (1.0-r[:,0])*6371., 'r-', linewidth=1.0)

#line_ref, = plt.plot(vsv_ref, (1.0-r_ref)*6371.0, 'b-', linewidth=1.0)

#line_TNA,    = plt.plot(vs_TNA, dep_TNA, 'm-', linewidth=0.5)
#line_TNA2,   = plt.plot(vs_TNA2, dep_TNA2, 'm.-', linewidth=0.5)
line_stw105, = plt.plot(vs_stw105, dep_stw105, 'k-', linewidth=1.5)

line_ref, = plt.plot(vs_ref, dep_ref, 'b-', linewidth=1.5)
line_FWEA18_1d_mean, = plt.plot(vs_FWEA18_1d_mean, dep_FWEA18_1d_mean, 'b-.', linewidth=1.0)
#line_prem,   = plt.plot(vs_prem, dep_prem, 'm-', linewidth=1.5)
#line_ak135,  = plt.plot(vs_ak135, dep_ak135, 'c-', linewidth=1.5)

#plt.legend([line_median, line_stw105, line_prem, line_ak135],
#    ['median', 'STW105', 'PREM', 'AK135F'])

plt.legend([line_mean, line_stw105, line_ref, line_FWEA18_1d_mean ],
    ['EARA2014_mean', 'STW105', 'FWEA18_REF', 'FWEA18_mean'])

#plt.legend([line_TNA2], ['TNA2 + 0.1'])

plt.xlim([4.5, 6.5])
plt.ylim([200, 800])
#plt.xlim([4.6, 5.4])
#plt.ylim([350, 480])

plt.grid()

plt.xlabel('Velocity (km/s)')
plt.ylabel('Depth (km)')

plt.gca().invert_yaxis()

#plt.show()
plt.savefig(out_fig)