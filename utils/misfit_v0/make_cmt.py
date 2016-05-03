#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit

import numpy as np
import matplotlib.pyplot as plt

# read command line args
misfit_file = "misfit/misfit.pkl"
cmt_file = "CMTSOLUTION.grid_cc"

#------
print "\n====== initialize\n"
misfit = Misfit()

print "\n====== load data\n"
misfit.load(filename=misfit_file)

# get best model
event = misfit.data['event']
src_perturb = misfit.data['src_perturb']

model_alpha = {
  'tau': 2.67503e-02,
   'mt': 1.57233e-01,
   't0': 1.20793e-01,
   'xs': 1.18584e+00,
  }

model = {}
for par in model_alpha:
  m0 = event[par]
  dm = src_perturb[par]
  alpha = model_alpha[par]
  model[par] = {
      'm0':m0, 
      'dm':dm,
      'alpha':alpha,
      'm1':m0+alpha*dm
      }

print model

#------ write out new CMTSOLUTION file
t0 = event['t0']
if 't0' in model: t0 = model['t0']['m1']
# modify origin time in header line to have centroid time 
header = event['header'].split()
header[1] = str(t0.year)
header[2] = str(t0.month)
header[3] = str(t0.day)
header[4] = str(t0.hour)
header[5] = str(t0.minute)
header[6] = str(t0.second + 1.0e-6*t0.microsecond)

tau = event['tau']
if 'tau' in model: tau = model['tau']['m1']
xs = event['xs']
if 'xs' in model: xs = model['xs']['m1']
mt = event['mt']
m0 = (0.5*np.sum(mt**2))**0.5
if 'mt' in model: mt = model['mt']['m1']
mt = mt - np.identity(3)*np.trace(mt)/3.0
mt = mt/(0.5*np.sum(mt**2))**0.5
mt *= m0

with open(cmt_file, 'w') as fp:
  fp.write('%s \n' % (' '.join(header)))
  fp.write('%-18s %s\n' % ('event name:', "grid_cc"))
  fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
  fp.write('%-18s %+15.8E\n' % ('tau(s):',   tau))
  fp.write('%-18s %+15.8E\n' % ('x(m):',     xs[0]))
  fp.write('%-18s %+15.8E\n' % ('y(m):',     xs[1]))
  fp.write('%-18s %+15.8E\n' % ('z(m):',     xs[2]))
  fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mt[0,0]))
  fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', mt[1,1]))
  fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mt[2,2]))
  fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mt[0,1]))
  fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mt[0,2]))
  fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', mt[1,2]))
