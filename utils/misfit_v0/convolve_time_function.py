#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convolve source time function
"""
import sys
import numpy as np
from obspy import read

#======
def stf_spectrum_gauss(f, tau):
  """ spectrum of the Gaussian STF of unit area:
      stf(t,tau) = 1/sqrt(PI)/tau * exp(-(t/tau)^2)
      F_stf = exp(- pi^2 * f^2 * tau^2)
  """
  return np.exp(-np.pi**2 * f**2 * tau**2)

#====== main
in_dir = str(sys.argv[1])
out_dir = str(sys.argv[2])
tau = float(sys.argv[3])

st = read(in_dir + '/*.sac')

for tr in st:
  n = tr.stats.npts
  dt = tr.stats.delta

  # pad tailing zeros to avoid signal wrap round contaminating the beginning
  nz = n + int(np.ceil(5*tau/dt)) # 5 is chosen based on full width of STF
  x = np.zeros(nz)
  x[0:n] = tr.data

  f = np.fft.rfftfreq(nz, d=dt)
  F_src = stf_spectrum_gauss(f, tau)
  x = np.fft.irfft(F_src*np.fft.rfft(x), nz)

  tr.data = x[0:n]

  tr.write(out_dir + '/' + tr.id, format="sac")

#END