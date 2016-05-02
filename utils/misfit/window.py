#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os 
import warnings
import re
#
import numpy as np
import scipy.signal as signal
#
import pickle
#
from obspy import UTCDateTime, read, Trace, geodetics
from obspy.taup import TauPyModel
from obspy.imaging.beachball import Beach
#
import pyproj
#
from lanczos_interp1 import lanczos_interp1
#
from matplotlib import colors, ticker, cm
from matplotlib.backends.backend_pdf import PdfPages 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#
from taper import *

#NOTE
# 1. spectrum relation between DFT and FT
#   x(n*dt): DFT[x]*dt ~ FT[x], IDFT[FT[x]]/dt ~ x

#====== utility
def is_equal(lst):
  return not lst or [lst[0]]*len(lst) == lst

def stf_gauss_spectrum(f, tau):
  """ 
  Spectrum of unit area gaussian function:
  s(t,tau) = 1/sqrt(PI)/tau * exp(-(t/tau)^2)
  FT(s) = exp(- pi^2 * f^2 * tau^2)

  f : 1d array of frequencies (Hz)
  tau : scalar, gauss width (sec)
  """
  return np.exp(-np.pi**2 * f**2 * tau**2)

#======
class Window:
  """
  time window for misifit/adjoint source

  Attributes
  ----------
  code
  starttime
  endtime
  taper
  taper_percentage
  filter
  filter_order
  filter_flo
  filter_fhi
  polarity_proj: (3,3) 

  quality{
    SNR, Amax_obs, Amax_noise, Amax_syn
  }

  cc{
    cc, dt_cc, CC0, CCmax, AR0, ARmax
  }

  weight

  Methods
  -------
  _shift_green_to_data
  measure_cc
  make_adjoint_source
  grid_search_cc

  """
#======================================================
  def __init__(self, starttime, endtime):
    self.status = 0
    self.history = ""

    self.label = "None"
    self.projection = "F"

    self.starttime = starttime
    self.endtime = endtime

    self.taper = 'cosine'
    self.taper_percentage = 0.1

    self.filter = 'butter'
    self.filter_order = 2
    self.filter_flo = 0.01
    self.filter_fhi = 0.1

    self.cc = {}
    self.quality = {}
    self.weight = 1.0

    # internal variables
    self._taper_func = None
    self._filter_a = None
    self._filter_b = None
    self._proj_matrix = None

#======================================================
  def _update_status(self, status, msg):
    self.status = status
    msg = "[{:s}] {:s}\n".format(UTCDateTime.now().isoformat(), msg)
    self.history += msg

#======================================================
  def _setup_window(self,
      sampling_rate,
      data_starttime, data_npts, back_azimuth,
      ):
    """ 
    setup window parameters: filter, taper

    Parameters
    ----------
    back_azimuth : scalar

    """
    #------ check window time range inside syn time
    data_endtime = data_starttime + data_npts/sampling_rate
    window_starttime = self.starttime
    window_endtime = self.endtime

    if window_endtime > data_endtime:
      msg = "window_endtime(%s) > data_endtime(%s)" % (
          window_endtime, data_endtime)
      raise Exception(msg)

    if window_starttime < data_starttime:
      msg = "window_starttime(%s) < data_starttime(%s)" % (
          window_endtime, data_endtime) 
      raise Exception(msg)

    #------ taper design
    win_b = window_starttime - data_starttime
    win_e = window_endtime - data_starttime
    win_len = window_endtime - window_starttime
    taper_width = win_len * min(taper_percentage, 0.5)
    win_c = [win_b, win_b+taper_width, win_e-taper_width, win_e]
    data_times = np.arange(data_npts)/sampling_rate
    if self.taper == "cosine":
      self._taper_func = cosine_taper(data_times, win_c)
    else:
      msg = "taper({:s}) not recognized.".format(self.taper)
      raise Exception(msg)

    #------ filter design
    nyq_freq = 0.5/sampling_rate
    freqlim = np.array([self.filter_flo, self.filter_fhi])/nyq_freq
    if self.filter == "butter":
      b, a = signal.butter(self.filter_order, freqlim, btype='band')
    else:
      msg = "filter({:s}) not recognized".format(self.filter)
      raise Exception(msg)
    self._filter_a = a
    self._filter_b = b

    #------ projection matrix
    comp = self.projection
    if comp == 'Z': # vertcal component
      cmpaz = 0.0 
      cmpdip = -90.0
    elif comp == 'R': # radial component
      cmpaz = (back_azimuth + 180.0)%360.0
      cmpdip = 0.0
    elif comp == 'T': # tangential component (TRZ: right-hand convention)
      cmpaz = (back_azimuth - 90.0)%360.0
      cmpdip = 0.0
    elif comp in ['H', 'I']: # Horizontal or Full  
      pass
    else:
      msg = "unrecognized projection {:s}".format(comp)
      raise Exception(msg)

    self._proj_matrix = np.zeros((3,3))
    if comp in ['Z', 'R', 'T']:
      sin_az = np.sin(np.deg2rad(cmpaz))
      cos_az = np.cos(np.deg2rad(cmpaz))
      sin_dip = np.sin(np.deg2rad(cmpdip))
      cos_dip = np.cos(np.deg2rad(cmpdip))
      n = np.array([
        [cos_dip * sin_az], # cos(E, v), v: projection vector
        [cos_dip * cos_az], # N
        [-sin_dip]] )     # Z
      self._proj_matrix = np.dot(n, n.transpose())
    elif comp == 'H': # horizontal
      self._proj_matrix[0,0] = 1.0 # E
      self._proj_matrix[1,1] = 1.0 # N
      self._proj_matrix[2,2] = 0.0 # Z
    elif comp == 'I': # Identity
      self._proj_matrix = np.identity(3)

#======================================================
  def _convolve_stf(self, sampling_rate, green, t0, gauss_width):
    """
    Make synthetic seismograms by convolving green's function with source time function

    Parameters
    ----------
    t0 : scalar, optional
    tau : scalar, optional
    """
    #------ source spectrum
    npts = green.shape[1]
    freq = np.fft.rfftfreq(npts, d=1.0/sampling_rate)
    src_spectrum = stf_gauss_spectrum(freq, gauss_width)

    #------ time shift
    ph_shift = np.exp(-2.0j*np.pi*freq*t0)

    #------ convolve STF
    return np.fft.irfft(src_spectrum*np.fft.rfft(green)*ph_shift, npts)

#======================================================
  def _process_data(self, 
      sampling_rate, data_starttime, data, green,
      t0, gauss_width, first_arrtime, back_azimuth,
      ):
    """ 
    Process data (filter, taper, projection, convlove stf)

    Parameters
    ----------
    data : ndarray (3,n)
    green : ndarray (3,n)
      green's function shifted/cut/pad to the same time samples as data
      with impulsive source at centroid_time
    t0 : scalar
      relative time shift of centroid_time
    gauss_width : scalar, > 0
      source time function parameter (gaussian)

    """
    #------ setup window parameters
    data_npts = data.shape[1]
    self._setup_window(sampling_rate, data_starttime, data_npts, back_azimuth)

    #------ filter obs, syn
    #NOTE: use lfilter (causal filter) to avoid contamination from the right
    # end of the signal, but with asymmetric response and 
    # peak shift ~ 1/4 min. period (e.g. 0.01-0.1Hz -> 2.5s peak shift)
    # , however the duration of the filter response is determined by the
    # max. period (e.g. 0.01-0.1Hz -> ~50s). So the time window chosen 
    # should not be affected by the relatively small peak shift.
    #-- F * d
    data_filt = signal.lfilter(self._filter_b, self._filter_a, data)
    #-- syn = stf * green
    syn = self._convolve_stf(sampling_rate, green, t0, gauss_width)
    #-- F * u
    syn_filt = signal.lfilter(self._filter_b, self._filter_a, syn)
    #DEBUG
    #for i in range(3):
    #  plt.subplot(311+i)
    #  plt.plot(syn_times, obs[i,:], 'k', syn_times, obs_filt[i,:], 'r')
    #plt.show()
    #-- noise: use signals 40s before first arrival time on data
    #FIXME: better choice of the time length before first arrival? 
    noise_idx1 = int(((first_arrtime - data_starttime) - 30.0)*sampling_rate)
    #t = syn_times[noise_idx]
    #b = t[0]
    #e = t[-1]
    #taper_width = (e-b) * 0.1
    #win_c = [b, b+taper_width, e-taper_width, e]
    #taper = cosine_taper(t, win_c)
    # F * noise
    noise_filt = data_filt[:,0:noise_idx1]

    #------ apply taper window and projection
    # w * F * d
    data_filt_win = np.dot(self._proj_matrix, data_filt) * self._taper_func
    # w * F * u
    syn_filt_win = np.dot(self._proj_matrix, syn_filt) * self._taper_func
    # noise (only projection)
    noise_filt_win = np.dot(self._proj_matrix, noise_filt)
    #DEBUG
    #diff = obs_ENZ_win - syn_ENZ_win
    #for i in range(3):
    #  plt.subplot(311+i)
    #  plt.plot(syn_times, obs_ENZ_win[i,:], 'k')
    #  plt.plot(syn_times, syn_ENZ_win[i,:], 'r')
    #  plt.plot(syn_times, diff[i,:], 'c')
    #plt.show()

    return data_filt_win, syn_filt_win, noise_filt_win

#======================================================
  def measure_cc(self, 
      sampling_rate, data_starttime, data, green,
      t0, gauss_width, first_arrtime, back_azimuth,
      plot=False,
      cc_delta=0.01, 
      ):
    """ 
    Measure data misfit based on cross-correlation value

    Parameters
    ----------
    data : ndarray (3,n)
    green : ndarray (3,n)
      green's function shifted/cut/pad to the same time samples as data
      with impulsive source at centroid_time
    t0 : scalar
      relative time shift of centroid_time
    gauss_width : scalar, > 0
      source time function parameter (gaussian)

    """
    data_filt_win, syn_filt_win, noise_filt_win = self._process_data(
        sampling_rate, data_starttime, data, green, 
        t0, gauss_width, first_arrtime, back_azimuth)

    #------ measure SNR (based on maximum amplitude)
    Amax_data = np.sqrt(np.max(np.sum(data_filt_win**2, axis=0)))
    Amax_syn = np.sqrt(np.max(np.sum(syn_filt_win**2, axis=0)))
    Amax_noise =  np.sqrt(np.max(np.sum(noise_filt_win**2, axis=0)))
    # empty data record
    if Amax_data == 0: 
      msg = "empty data trace."
      raise Exception(msg)
    # could occure when the data begin time is too close to the first arrival
    if Amax_noise == 0: 
      msg = "empty noise trace."
      raise Exception(msg)
    # SNR
    snr = 20.0*np.log10(Amax_data/Amax_noise)
 
    #------ measure CC time shift
    data_norm = np.sqrt(np.sum(data_filt_win**2))
    syn_norm = np.sqrt(np.sum(syn_filt_win**2))
    # window normalization factor (without *dt)
    Nw = data_norm * syn_norm
    cc[:] = 0.0
    # NOTE the order (obs,syn) is important. The positive time on 
    # CC means shifting syn along the positive time direction in order to match
    # the observed obs, and vice verser.
    # [-(nt-1), nt) * dt
    for i in range(3):
      cc += signal.fftconvolve(
          data_filt_win[i,:], syn_filt_win[i,::-1], 'full')
    cc /= Nw
    #DEBUG
    #print window_id
    #print cc[syn_nt-2] - np.sum(obs_ENZ_win * syn_ENZ_win)
    #print cc[syn_nt-1] - np.sum(obs_ENZ_win * syn_ENZ_win)
    #print cc[syn_nt] - np.sum(obs_ENZ_win * syn_ENZ_win)
    #print cc[syn_nt+1] - np.sum(obs_ENZ_win * syn_ENZ_win)
    #print cc[syn_nt+2] - np.sum(obs_ENZ_win * syn_ENZ_win)
    #-- zero-lag cc coeff.
    data_npts = data.shape[1]
    CC0 = cc[data_npts-1] #the n-th point corresponds to zero lag time 
    AR0 = CC0 * syn_norm / data_norm # amplitude ratio syn/obs 
    #DEBUG
    #print CC0 - np.sum(obs_ENZ_win * syn_ENZ_win)/obs_norm/syn_norm
    #-- interpolate cc to finer time samples
    win_len = self.endtime - self.starttime
    CC_shift_range = win_len/2.0 #TODO: more reasonable choice?
    ncc = int(CC_shift_range / cc_delta)
    cc_times = np.arange(-ncc,ncc+1) * cc_delta
    delta = 1.0/sampling_rate
    if delta < cc_delta:
      msg = "delta(%f) < cc_time_step(%f)" % (delta, cc_delta)
      warnings.warn(msg)
    ti = cc_times + (data_npts-1)*delta  # -(npts-1)*dt: begin time in cc
    cci = lanczos_interp1(cc, delta, ti, na=20)
    # time shift at the maximum correlation
    imax = np.argmax(cci)
    CC_time_shift = cc_times[imax]
    CCmax = cci[imax]
    ARmax = CCmax * syn_norm / obs_norm # amplitude ratio: syn/obs

    # adjoint source: dchiw_du (misfit functional: zero-lag cc coef.)
    # dchiw_du = conj(F * [S]) * w * [ w * F * d - A * w * F * S * g] / N, 
    # , where A = CC0(un-normalized) / norm(u)**2, N = norm(d)*norm(u)
    Aw = CC0 * data_norm / syn_norm # window amplitude raito

    #------ record results
    quality_dict = {
        'Amax_data': Amax_data, 'Amax_syn': Amax_syn, 
        'Amax_noise': Amax_noise, 'SNR': snr}

    cc_dict = {
        'time': cc_times, 'cc': cci,
        'time_shift': CC_time_shift,
        'CC0': CC0, 'CCmax': CCmax,
        'AR0': AR0, 'ARmax': ARmax,
        'Nw':Nw, 'Aw':Aw }

    self.quality.update(quality_dict)
    self.cc.update(cc_dict)

    msg = "measure_cc OK."
    self._update_status(0, msg)

#======================================================
  def make_adjoint_source(self, 
      adj_type='dchi_dg',
      out_dir='adj',
      band_code='MX'):
    """
    Output adjoint source

    NOTE
    ----
    dchi_dg: use tau=0 in forward/adjoint simulation

    dchi_du: use real tau in forward/adjoint simulation

    """
    orientation_codes = ['E', 'N', 'Z']
    event = self.data['event']
    sampling_rate = self.data['sampling_rate']
    delta = 1.0/sampling_rate

    tr = Trace()
    station_dict = self.data['station']
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue
      # adjoint source
      if adj_type == 'dchi_du':
        adj = station['dchi_du']
      elif adj_type == 'dchi_dg':
        adj = station['dchi_dg']
      else:
        raise Exception('unknown adj_type: %s (dchi_du or dchi_dg) ' \
            % (adj_type))

      # time samples
      waveform = station['waveform']
      obs_starttime = waveform['obs_starttime']
      grf_t0 = waveform['grf_t0']
      simulation_delta = waveform['simulation_delta']
      simulation_npts = waveform['simulation_npts']
      simulation_times = np.arange(simulation_npts)*simulation_delta + grf_t0

      # loop ENZ
      time_shift = event['t0'] - obs_starttime
      for i in range(3):
        # interpolate back to synthetic time samples (required by SEM)
        adj1 = lanczos_interp1(adj[i,:], delta, 
            simulation_times + time_shift, na=20)
        tr.data = adj1
        tr.stats.starttime = event['t0'] + syn_t0
        tr.stats.delta = syn_delta
        out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
            out_dir, station_id, band_code,
            orientation_codes[i])
        # sac format
        tr.write(out_file + '.adj.sac', 'sac')
        # ascii format (needed by SEM)
        # time is relative to event origin time: t0
        with open(out_file+'.adj','w') as fp:
          for j in range(simulation_npts):
            fp.write("{:16.9e}  {:16.9e}\n".format(
              simulation_times[j], adj1[j]))

    #endfor station_id in station_dict:
  #enddef make_adjoint_source

#
#======================================================
#

  def read_srcfrechet(self, filename=None, update=False):
    """ Read in source derivative of misfit function
        Dchi/Dxs, Dchi/Dmt
    """
    with open(filename, 'r') as f:
      lines = [ x for x in f.readlines() if not(x.startswith('#')) ]

    lines = [x.split() for x in lines]

    t0  = float(lines[0][0]);  dchi_dt0  = float(lines[0][1])
    tau = float(lines[1][0]);  dchi_dtau = float(lines[1][1])
    x   = float(lines[2][0]);  dchi_dx   = float(lines[2][1])
    y   = float(lines[3][0]);  dchi_dy   = float(lines[3][1])
    z   = float(lines[4][0]);  dchi_dz   = float(lines[4][1])
    mxx = float(lines[5][0]);  dchi_dmxx = float(lines[5][1])
    myy = float(lines[6][0]);  dchi_dmyy = float(lines[6][1])
    mzz = float(lines[7][0]);  dchi_dmzz = float(lines[7][1])
    mxy = float(lines[8][0]);  dchi_dmxy = float(lines[8][1])
    mxz = float(lines[9][0]);  dchi_dmxz = float(lines[9][1])
    myz = float(lines[10][0]); dchi_dmyz = float(lines[10][1])

    dchi_dxs = np.array([dchi_dx, dchi_dy, dchi_dz])

    dchi_dmt = np.zeros((3,3))
    dchi_dmt[0,0] = dchi_dmxx
    dchi_dmt[1,1] = dchi_dmyy
    dchi_dmt[2,2] = dchi_dmzz
    dchi_dmt[0,1] = dchi_dmxy
    dchi_dmt[1,0] = dchi_dmxy
    dchi_dmt[0,2] = dchi_dmxz
    dchi_dmt[2,0] = dchi_dmxz
    dchi_dmt[1,2] = dchi_dmyz
    dchi_dmt[2,1] = dchi_dmyz

    # check if the same as event info
    data = self.data
    event = data['event']
    #...

    # record 
    src_frechet = {
        't0':dchi_dt0,
        'tau':dchi_dtau,
        'xs':dchi_dxs,
        'mt':dchi_dmt,
        'stat': {'code':0, 'msg':"created on "+UTCDateTime.now().isoformat()}
        }

    if 'src_frechet' not in data:
      data['src_frechet'] = src_frechet
    elif update:
      data['src_frechet'].update(src_frechet)
      data['src_frechet']['stat']['code'] = 1
      data['src_frechet']['stat']['msg'] = "updated on "+UTCDateTime.now().isoformat()
    else:
      raise Exception('src_frechet already set, not updated.')

#
#======================================================
#

  def make_cmt_dxs(self, out_file="CMTSOLUTION.dxs", norm=2500.0):
    """ Calculate derivative for source location along one direction
    """
    norm = float(norm)
    if norm <= 0.0:
      raise Exception("norm(dxs) must be larger than 0.")

    # initialize pyproj objects
    geod = pyproj.Geod(ellps='WGS84')
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    # get source parameters
    event = self.data['event']
    tau = event['tau']
    xs = event['xs']
    mt = event['mt']

    # get perturbed source location
    if 'src_frechet' not in self.data:
      raise Exception('src_frechet is not set.')
    src_frechet = self.data['src_frechet']
    dxs = src_frechet['xs']
    # normalize dxs
    dxs = dxs/(np.sum(dxs**2))**0.5
    # apply given norm 
    dxs *= norm
    # get new src location
    xs1 = xs + dxs
    lon, lat, alt = pyproj.transform(ecef, lla, xs1[0], xs1[1], xs1[2])
    depth = -alt
    if depth < 0.0:
      raise Exception("new src depth %f < 0.0" % depth)

    # record dxs
    if 'src_perturb' not in self.data:
      self.data['src_perturb'] = {}
    self.data['src_perturb']['xs'] = dxs

    # write out new CMTSOLUTION file
    with open(out_file, 'w') as fp:
      fp.write('%s\n' % event['header'])
      fp.write('%-18s %s_dxs\n' % ('event name:',event['id']))
      fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
      fp.write('%-18s %+15.8E\n' % ('tau(s):',   0.0))
      fp.write('%-18s %+15.8E\n' % ('x(m):',     xs1[0]))
      fp.write('%-18s %+15.8E\n' % ('y(m):',     xs1[1]))
      fp.write('%-18s %+15.8E\n' % ('z(m):',     xs1[2]))
      fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mt[0,0]))
      fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', mt[1,1]))
      fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mt[2,2]))
      fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mt[0,1]))
      fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mt[0,2]))
      fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', mt[1,2]))

#
#======================================================
#

  def waveform_der_dxs(self,
      syn_dir='output_dxs',
      syn_band_code='MX',
      syn_suffix='.sem.sac',
      sac_dir=None):
    """ Calculate derivative for source location along one direction
    NOTE:
      1) use finite difference to get waveform derivative
      2) dxs: length 3 vector (unit: meter)
      3) use green's function as input (i.e. set tau to zero in simulation)
    """
    syn_orientation_codes = ['E', 'N', 'Z']

    event = self.data['event']
    tau = event['tau']
    t0 = event['t0']

    # src_perturb
    dxs = self.data['src_perturb']['xs']

    station_dict = self.data['station']
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue

      #------ time samples
      waveform = station['waveform']
      data_starttime = time_sample['starttime']
      dt = time_sample['delta']
      nt = time_sample['nt']
      nl = time_sample['nl'] # npts of left padding
      nr = time_sample['nr'] # npts of right padding
      t = np.arange(nt) * dt + (starttime - t0) #referred to t0

      #------ get file paths of syn seismograms
      syn_files = [ '{:s}/{:s}.{:2s}{:1s}{:s}'.format(
        syn_dir, station_id, syn_band_code, x, syn_suffix)
        for x in syn_orientation_codes ]

      #------ read in obs, syn seismograms
      syn_st  = read(syn_files[0])
      syn_st += read(syn_files[1])
      syn_st += read(syn_files[2])

      #------ check the same time samples as original syn
      if not is_equal( [ (tr.stats.starttime, tr.stats.delta, tr.stats.npts) \
          for tr in syn_st ] ):
        raise Exception('%s: not equal time samples in'\
            ' synthetic seismograms.' % (station_id))
      tr = syn_st[0]
      if tr.stats.delta != dt:
        raise Exception("%s: not the same dt for diff-srcloc!" % (station_id))
      if (tr.stats.starttime - nl*dt) != starttime:
        raise Exception("%s: not the same origin time for diff-srcloc!" % (station_id))
      if tr.stats.npts != (nt-nl-nr):
        raise Exception("%s: not the same npts for diff-srcloc!" % (station_id))

      #------ read syn seismograms from perturbed source location
      syn_ENZ = np.zeros((3, nt))
      for i in range(3):
        syn_ENZ[i,nl:(nl+nt)] = syn_st[i].data

      # differential green's function 
      grf = waveform['grf']
      dg = syn_ENZ - grf

      ## diff synthetics
      ## convlove source time function
      #syn_freq = np.fft.rfftfreq(nt, d=dt)
      #F_src = stf_gauss_spectrum(syn_freq, tau)
      #du = np.fft.irfft(F_src * np.fft.rfft(dg), nt)
      ##zero records before origin time (wrap around from the end)
      #idx = t < -5.0*tau
      #du[:,idx] = 0.0

      ## diff Chi
      #dchi = np.sum(station['dchi_dg'] * dg) * dt

      #------ record derivatives
      if 'waveform_der' not in station:
        station['waveform_der'] = {}
      station['waveform_der']['xs'] = {
          'dm':dxs, 'dg':dg }
      #'dm':dxs, 'dg':dg, 'du':du, 'dchi':dchi }

      # DEBUG
      #print dchi
      #for i in range(3):
      #  plt.subplot(311+i)
      #  #plt.plot(t,grf0[i,:],'k', t,syn_ENZ[i,:],'r', t,dg[i,:], 'b')
      #  plt.plot(t, du[i,:], 'k')
      #plt.show()
      if sac_dir:
        for i in range(3):
          tr.data = du[i,:]
          tr.stats.starttime = starttime
          tr.stats.delta = dt
          tr.stats.npts = nt
          out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
              sac_dir, station_id, syn_band_code,
              syn_orientation_codes[i])
          tr.write(out_file, 'sac')

#
#======================================================
#

  def make_cmt_dmt(self,
      out_file="CMTSOLUTION.dmt",
      fix_M0=True, zerotrace=True, ratio_M0=0.01):
    """ Calculate derivative for source location along one direction
      fix_M0: project dmt orthogonal mt to keep seismic moment M0 = sqrt(0.5*m:m) fixed
      zerotrace: zero tr(dmt)
      ratio_M0: magnitude of dmt is set to the ratio of M0
    """
    # get source parameters
    event = self.data['event']
    tau = event['tau']
    xs = event['xs']
    mt = event['mt']
    
    # check parameters
    if ratio_M0 <= 0.0:
      error_str = "ratio_M0(%f) must > 0" % ratio_M0
      raise ValueError(error_str)

    # get perturbed moment tensor
    if 'src_frechet' not in self.data:
      raise Exception('src_frechet not set.')
    src_frechet = self.data['src_frechet']
    # set dmt parallel to dchi_dmt
    dmt = src_frechet['mt']
    # project dmt perpendicular to M0 change direction
    if fix_M0:
      dmt = dmt - mt*np.sum(dmt*mt)/np.sum(mt**2)
    if zerotrace:
      dmt = dmt - np.identity(3)*np.trace(dmt)/3.0
    # normalize dmt to have unit seismic moment
    dmt = dmt/(0.5*np.sum(dmt**2))**0.5
    # use ratio_M0 as the magnitude of dmt
    m0 = (0.5*np.sum(mt**2))**0.5
    dmt *= ratio_M0 * m0

    # record dmt
    if 'src_perturb' not in self.data:
      self.data['src_perturb'] = {}
    self.data['src_perturb']['mt'] = dmt

    # write out new CMTSOLUTION file
    with open(out_file, 'w') as fp:
      fp.write('%s\n' % event['header'])
      fp.write('%-18s %s_dmt\n' % ('event name:',event['id']))
      fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
      fp.write('%-18s %+15.8E\n' % ('tau(s):',   0.0))
      fp.write('%-18s %+15.8E\n' % ('x(m):',     xs[0]))
      fp.write('%-18s %+15.8E\n' % ('y(m):',     xs[1]))
      fp.write('%-18s %+15.8E\n' % ('z(m):',     xs[2]))
      fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', dmt[0,0]))
      fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', dmt[1,1]))
      fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', dmt[2,2]))
      fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', dmt[0,1]))
      fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', dmt[0,2]))
      fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', dmt[1,2]))

#
#======================================================
#

  def waveform_der_dmt(self,
      syn_dir='output_dmt',
      syn_band_code='MX',
      syn_suffix='.sem.sac',
      sac_dir=None):
    """ Calculate derivative for moment tensor along a given direction
    NOTE:
      1) dmt: 3 by 3 symetric matrix (unit: N*m)
      2) use green's function as input (i.e. set tau to zero in simulation)
    """
    syn_orientation_codes = ['E', 'N', 'Z']
    # event
    event = self.data['event']
    tau = event['tau']
    t0 = event['t0']

    # src_perturb
    dmt = self.data['src_perturb']['mt']

    station_dict = self.data['station']
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue

      #------ time samples
      waveform = station['waveform']
      time_sample = waveform['time_sample']
      starttime = time_sample['starttime']
      dt = time_sample['delta']
      nt = time_sample['nt']
      nl = time_sample['nl'] # npts of left padding
      nr = time_sample['nr'] # npts of right padding
      t = np.arange(nt) * dt + (starttime - t0) #referred to t0

      #------ get file paths of syn seismograms
      syn_files = [ '{:s}/{:s}.{:2s}{:1s}{:s}'.format(
        syn_dir, station_id, syn_band_code, x, syn_suffix)
        for x in syn_orientation_codes ]

      #------ read in obs, syn seismograms
      syn_st  = read(syn_files[0])
      syn_st += read(syn_files[1])
      syn_st += read(syn_files[2])

      #------ check the same time samples as original syn
      if not is_equal( [ (tr.stats.starttime, tr.stats.delta, tr.stats.npts) \
          for tr in syn_st ] ):
        raise Exception('%s: not equal time samples in'\
            ' synthetic seismograms.' % (station_id))
      tr = syn_st[0]
      if tr.stats.delta != dt:
        raise Exception("%s: not the same dt for diff-srcloc!" % (station_id))
      if (tr.stats.starttime - nl*dt) != starttime:
        raise Exception("%s: not the same origin time for diff-srcloc!" % (station_id))
      if tr.stats.npts != (nt-nl-nr):
        raise Exception("%s: not the same npts for diff-srcloc!" % (station_id))

      #------ read syn seismograms from perturbed source location
      dg = np.zeros((3, nt))
      for i in range(3):
        dg[i,nl:(nl+nt)] = syn_st[i].data

      ##source spectrum (moment-rate function)
      #syn_freq = np.fft.rfftfreq(nt, d=dt)
      #F_src = stf_gauss_spectrum(syn_freq, tau)
      #du = np.fft.irfft(F_src * np.fft.rfft(dg), nt)
      ##zero records before origin time (wrap around from the end)
      #idx = t < -5.0*tau
      #du[:,idx] = 0.0

      ## diff Chi
      #dchi = np.sum(station['dchi_dg'] * dg)

      #------ record derivatives
      if 'waveform_der' not in station:
        station['waveform_der'] = {}
      station['waveform_der']['mt'] = {
          'dm':np.array(dmt), 'dg':dg }
      #'dm':np.array(dmt), 'dg':dg, 'du':du, 'dchi':dchi }

      # DEBUG
      #print dchi
      #for i in range(3):
      #  plt.subplot(311+i)
      #  plt.plot(t, du[i,:], 'k')
      #plt.show()
      if sac_dir:
        for i in range(3):
          tr.data = du[i,:]
          tr.stats.starttime = starttime
          tr.stats.delta = dt
          tr.stats.npts = nt
          out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
              sac_dir, station_id, syn_band_code,
              syn_orientation_codes[i])
          tr.write(out_file, 'sac')

#
#======================================================
#

  def measure_hessian_src(self, update=False):
    """ calculate hessian matrix for source parameters (dchi_du)
        chi: misfit functional (normalized zero-lag correlation coef.)
        u: synthetic waveform
    """
    event = self.data['event']
    src_param = ('dt0','dtau','dxs','dmt')
    n_srcparam = len(src_param)

    #------ loop each station
    station_dict = self.data['station']
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue

      # waveform
      waveform = station['waveform']
      time_sample = waveform['time_sample']
      syn_starttime = time_sample['starttime']
      syn_delta = time_sample['delta']
      syn_nyq = 0.5/syn_delta
      syn_nt = time_sample['nt']
      syn_nl = time_sample['nl']
      syn_nr = time_sample['nr']
      syn_times = syn_delta * np.arange(syn_nt)
      # seismograms 
      obs = waveform['obs']
      grf = waveform['grf']

      # source spectrum (moment-rate function)
      syn_freq = np.fft.rfftfreq(syn_nt, d=syn_delta)
      F_src = stf_gauss_spectrum(syn_freq, event['tau'])

      # waveform derivatives
      waveform_der = station['waveform_der']

      #------ loop each window
      window_dict = station['window']
      for window_id in window_dict:
        # window parameters
        window = window_dict[window_id]
        # skip bad windows
        if window['stat']['code'] < 1:
          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
          continue

        #------ window parameters 
        # filter
        filter_dict = window['filter']
        filter_a = filter_dict['a']
        filter_b = filter_dict['b']
        # taper
        win_func = window['taper']['win']
        # polarity projection 
        proj_matrix = window['polarity']['proj_matrix']

        #------ filter obs, syn
        # F * d
        obs_filt = signal.lfilter(filter_b, filter_a, obs)
        # F * u (u = S*grf)
        syn_filt = signal.lfilter(filter_b, filter_a, grf)
        syn_filt = np.fft.irfft(F_src*np.fft.rfft(syn_filt), syn_nt)
        # apply window taper and polarity projection
        # obs = w * F * d
        wFd = np.dot(proj_matrix, obs_filt) * win_func 
        # syn = w * F * u (u = S*grf)
        wFu = np.dot(proj_matrix, syn_filt) * win_func 
        # norm
        norm_wFd = np.sqrt(np.sum(wFd**2))
        norm_wFu = np.sqrt(np.sum(wFu**2))
        # window normalization factor
        Nw = norm_wFd * norm_wFu
        # window amplitude raito
        Aw = np.sum(wFd * wFu) / norm_wFu**2

        #DEBUG
        #print "Nw: %e %e" % (Nw, window['cc']['Nw'])
        #print "Aw: %e %e" % (Aw, window['cc']['Aw'])

        #------ filter differential seismograms (w * F * du_dm)
        wFdu = {}
        for param in src_param:
          du = waveform_der[param]['du']
          Fdu = signal.lfilter(filter_b, filter_a, du)
          wFdu[param] = np.dot(proj_matrix, Fdu) * win_func 

        #------ hessian src
        # chi: zero-lag correlation coef. between wFu and wFd
        # hessian: ddchi_dmdm
        hessian_src = {}
        for i in range(n_srcparam):
          for j in range(i, n_srcparam):
            par1 = src_param[i]
            par2 = src_param[j]
            wFdu1 = wFdu[par1]
            wFdu2 = wFdu[par2]
            wFdu1_wFdu2 = np.sum(wFdu1 * wFdu2)
            wFu_wFdu1 = np.sum(wFu * wFdu1)
            wFu_wFdu2 = np.sum(wFu * wFdu2)
            wFd_wFdu1 = np.sum(wFd * wFdu1)
            wFd_wFdu2 = np.sum(wFd * wFdu2)
            key12 = (par1, par2)
            hessian_src[key12] = ( \
                - Aw * wFdu1_wFdu2 \
                + ( 3.0 * Aw * wFu_wFdu1 * wFu_wFdu2 \
                             - wFu_wFdu1 * wFd_wFdu2 \
                             - wFu_wFdu2 * wFd_wFdu1 \
                  ) / norm_wFu**2 
                ) / Nw

        #------ record results
        if 'hessian_src' not in window:
          window['hessian_src'] = hessian_src
          window['stat'] = {'code': 2,
              'msg': "add hessian_src on "+UTCDateTime.now().isoformat()}
        elif update:
          window['hessian_src'].update(hessian_src)
          window['stat'] = {'code': 2,
              'msg': "update hessian_src on "+UTCDateTime.now().isoformat()}
        else:
          warnings.warn("hessian_src already set, nothing changed")
      # end for window_id in windows:
    # endfor station_id in station_dict:
  #enddef measure_windows_for_one_station(self,

#
#======================================================
#

# def update_source(self):
#   """ Update source parameters based on waveform derivatives and hessian
#   """
#   event = self.data['event']
#   src_param = ('dt0','dtau','dxs','dmt')
#   n_srcparam = len(src_param)

#   dchi_dm = np.zeros(n_srcparam)
#   hessian = np.zeros([n_srcparam,n_srcparam])

#   #------ get dchi_dm and Hessian 
#   #-- loop each station
#   station_dict = self.data['station']
#   for station_id in station_dict:
#     station = station_dict[station_id]
#     # skip rejected statations
#     if station['stat']['code'] < 0:
#       continue

#     # dchi_dm
#     for i in range(n_srcparam):
#       key = src_param[i]
#       dchi_dm[i] += station['waveform_der'][key]['dchi']

#     #-- loop each window
#     window_dict = station['window']
#     for window_id in window_dict:
#       # window parameters
#       window = window_dict[window_id]
#       # skip bad windows
#       if window['stat']['code'] < 1:
#         warnings.warn("Window %s not measured for adj, SKIP" % window_id)
#         continue
#       # 
#       weight = window['weight']
#       hessian_win = window['hessian_src']
#       for i in range(n_srcparam):
#         for j in range(i, n_srcparam):
#           par1 = src_param[i]
#           par2 = src_param[j]
#           key = (par1,par2)
#           hessian[i,j] += weight * hessian_win[key]
#     #end for window_id in windows:

#   #end for station_id in station_dict:

#   for i in range(n_srcparam):
#     for j in range(i+1, n_srcparam):
#         hessian[j,i] = hessian[i,j]

#   print "dchi_dm:"
#   print dchi_dm 

#   print "hessian:"
#   print hessian

#   print "====== 0:4:"
#   w, v = np.linalg.eigh(hessian, UPLO='U')
#   print w
#   print v
#   x, residual, rank, sigval = np.linalg.lstsq(hessian, -dchi_dm)
#   print " inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print "dt0: \n", x[0]
#   print "dtau:\n", x[1]
#   print "dxs: \n", x[2]*self.data['src_perturb']['xs'] 
#   print "dmt: \n", x[3]*self.data['src_perturb']['mt'] 

#   print "====== only 0:3"
#   h3 = hessian[0:3,0:3]
#   v3 = dchi_dm[0:3]
#   w, v = np.linalg.eigh(h3, UPLO='U')
#   print w
#   print v
#   x, residual, rank, sigval = np.linalg.lstsq(h3, -v3)
#   print "inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print "dt0: \n", x[0]
#   print "dtau:\n", x[1]
#   print "dxs: \n", x[2]*self.data['src_perturb']['xs'] 
#   #print "dmt: \n", x[3]*self.data['src_perturb']['dmt'] 

#   print "====== only 0:2"
#   h3 = hessian[0:2,0:2]
#   v3 = dchi_dm[0:2]
#   w, v = np.linalg.eigh(h3, UPLO='U')
#   print w
#   print v
#   x, residual, rank, sigval = np.linalg.lstsq(h3, -v3)
#   print "inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print "dt0: \n", x[0]
#   print "dtau:\n", x[1]

#   print "====== only 0,2"
#   idx = [0,2]
#   hess = hessian[idx,:][:,idx]
#   kernel = dchi_dm[idx]
#   print hess
#   eigs, eigv = np.linalg.eigh(hess, UPLO='U')
#   print eigs
#   print eigv
#   x, residual, rank, sigval = np.linalg.lstsq(hess, -kernel)
#   print "inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print "dt0: \n", x[0]
#   print "dxs:\n", x[1]

# #enddef measure_windows_for_one_station(self,

#
#======================================================
#

  def cc_perturbed_seismogram(self,
      dm={'t0':None, 'xs':None},
      plot=False
      ):
    """ calculate normalized zero-lag cc for perturbed seismograms from linear combination of waveform derivatives

    dm={<model_name>:<model_vector>, ...}
    model_name: ['dt0', 'dxs']
    model_vector: ndarray of the same length

    return cc_sum, weight_sum
    """
    # check model vectors in dm have the same length
    if (not dm) or len(dm) == 0:
      error_str = "dm must not be empty"
      raise Exception(error_str)
    model_num = len(dm)

    vector_size = []
    for model_name in dm:
      if dm[model_name].ndim != 1:
        error_str = "dm[%s] must be vector" % model_name
        raise Exception(error_str)
      vector_size.append(np.size(dm[model_name]))
    if not is_equal(vector_size) or vector_size[0] < 1:
      error_str = "vectors in dm must have the same non-zero length"
      raise Exception(error_str)
    vector_size = vector_size[0]

    # check parameters
    event = self.data['event']
    if 'tau' in dm:
      tau = dm['tau'] + event['tau']
      if any(tau <= 0):
        error_str = "dm['dtau'] has invalid values (event['tau']=%f)!" \
            % event['tau']
        raise Exception(error_str)

    #------ loop each station
    # sum of weighted normalized zero-lag cc at each model grid
    wcc_sum = np.zeros(vector_size)
    # sum of all windows' weighting
    weight_sum = 0.0

    station_dict = self.data['station']
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue

      # check if model parameter included in waveform_der
      waveform_der = station['waveform_der']
      for model_name in dm:
        if (model_name not in ['t0', 'tau']) and \
           (model_name not in waveform_der):
          error_str = "%s not in waveform_der of %s" % (model_name, station_id)
          raise Exception(error_str)

      #---- get seismograms: obs,grf 
      waveform = station['waveform']
      obs = waveform['obs']
      grf = waveform['grf']
      # time samples
      time_sample = waveform['time_sample']
      syn_starttime = time_sample['starttime']
      syn_delta = time_sample['delta']
      syn_nt = time_sample['nt']
      syn_nl = time_sample['nl']
      syn_nr = time_sample['nr']
      syn_freq = np.fft.rfftfreq(syn_nt, d=syn_delta)

      #---- measure misfit
      window_dict = station['window']
      for window_id in window_dict:
        window = window_dict[window_id]
        # skip bad windows
        if window['stat']['code'] < 1:
          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
          continue
        # window weight
        weight = window['weight']
        weight_sum += weight
        # filter
        filter_dict = window['filter']
        filter_a = filter_dict['a']
        filter_b = filter_dict['b']
        # taper
        win_func = window['taper']['win']
        win_starttime = window['taper']['starttime']
        win_endtime = window['taper']['endtime']
        # polarity projection 
        proj_matrix = window['polarity']['proj_matrix']
        #-- filter,project,taper obs
        # F * d
        obs_filt = signal.lfilter(filter_b, filter_a, obs)
        # w * p * F * d (window,project,filter)
        wpFd = np.dot(proj_matrix, obs_filt) * win_func 
        norm_wpFd = np.sqrt(np.sum(wpFd**2))
        #-- filter,project grf 
        # F * g
        grf_filt = signal.lfilter(filter_b, filter_a, grf)
        # p * F * g
        pFg = np.dot(proj_matrix, grf_filt)
        if plot:
          F_src = stf_gauss_spectrum(syn_freq, event['tau'])
          # S * F * g
          syn_filt = np.fft.irfft(F_src*np.fft.rfft(grf_filt), syn_nt)
          # w * p * S * F * g
          wpFu = np.dot(proj_matrix, syn_filt) * win_func
        #-- filter,project dg: pFdg
        pFdg = {}
        for model_name in dm:
          # exclude source time function
          if model_name not in ['t0', 'tau']:
            dg = waveform_der[model_name]['dg']
            dg_filt = signal.lfilter(filter_b, filter_a, dg)
            pFdg[model_name] = np.dot(proj_matrix, dg_filt)
        #-- misfit function: zero-lag cc
        for idx_model in range(vector_size):
          # perturbed grf: pFg1
          pFg1 = np.zeros((3,syn_nt))
          pFg1 += pFg
          for model_name in dm:
            # exclude source time function
            if model_name not in ['t0', 'tau']:
              pFg1 += dm[model_name][idx_model] * pFdg[model_name]
          # perturbed source time function
          dt0 = 0.0
          if 't0' in dm:
            dt0 = dm['t0'][idx_model]
          dtau = 0.0
          if 'tau' in dm:
            dtau = dm['tau'][idx_model]
          F_src = stf_gauss_spectrum(syn_freq, event['tau']+dtau)
          # perturbed syn: w * S * p * F * g1
          phase_shift = np.exp(-2.0j * np.pi * syn_freq * dt0)
          wpFu1 = np.fft.irfft(phase_shift*F_src*np.fft.rfft(pFg1), syn_nt) \
              * win_func
          norm_wpFu1 = np.sqrt(np.sum(wpFu1**2))
          Nw = norm_wpFd * norm_wpFu1
          #normalized cc between obs and perturbed syn
          cc_wpFd_wpFu1 = np.sum(wpFd*wpFu1) / Nw
          # weighted cc
          wcc_sum[idx_model] += weight * cc_wpFd_wpFu1
          #DEBUG
          if plot:
            syn_npts = syn_nt - syn_nl - syn_nr
            syn_orientation_codes = ['E', 'N', 'Z']
            syn_times = np.arange(syn_nt) * syn_delta
            Amax_obs = np.max(np.sum(wpFd**2, axis=0))**0.5
            Amax_syn = np.max(np.sum(wpFu**2, axis=0))**0.5
            Amax_syn1 = np.max(np.sum(wpFu1**2, axis=0))**0.5
            win_b = win_starttime - syn_starttime
            win_e = win_endtime - syn_starttime
            for i in range(3):
              plt.subplot(311+i)

              if i == 0:
                title_str = "%s.%s " % (station_id, window_id)
                for model_name in dm:
                  title_str += "%s:%.2f " \
                      % (model_name, dm[model_name][idx_model])
                title_str += "NCCwin:%.2f" % (cc_wpFd_wpFu1)
                plt.title(title_str)

              idx_plt = range(syn_nl,(syn_nl+syn_npts))
              plt.plot(syn_times[idx_plt], obs_filt[i,idx_plt]/Amax_obs, 
                  'k', linewidth=0.2)
              plt.plot(syn_times[idx_plt], syn_filt[i,idx_plt]/Amax_syn,
                  'b', linewidth=0.2)

              idx_plt = (win_b <= syn_times) & (syn_times <= win_e)
              plt.plot(syn_times[idx_plt], wpFd[i,idx_plt]/Amax_obs,
                  'k', linewidth=1.0)
              plt.plot(syn_times[idx_plt], wpFu[i,idx_plt]/Amax_syn,
                  'b', linewidth=1.0)
              plt.plot(syn_times[idx_plt], wpFu1[i,idx_plt]/Amax_syn1,
                  'r', linewidth=1.0)

              plt.xlim((syn_times[syn_nl], syn_times[syn_nl+syn_npts-1]))
              plt.ylim((-1.5, 1.5))
              plt.ylabel(syn_orientation_codes[i])
            plt.show()

      #end for window_id in window_dict:
    #end for station_id in station_dict:

    return wcc_sum, weight_sum

#
#======================================================
#

  def grid_cc_perturbed_seismogram(self,
      dm = {
        't0': np.linspace(-10,0,11), 
        'xs': np.linspace(-5,5,11) 
        },
      axes=[ ('t0','xs'), ('t0',), ('xs',) ],
      outfig="grid_cc.pdf",
      plot_seism=False,
      cmt_file="CMTSOLUTION.grid_cc"
      ):
    """ calculate misfit over 2D model grids based on perturbed seismograms 
    """
    # grid parameters
    model_num = len(dm)
    model_name = [x for x in dm]

    # check parameters 
    for xy in axes:
      if len(xy) < 1 or len(xy) > 2:
        error_str = "axes should have one or two elements"
        raise Exception(error_str)
      for axis in xy:
        if axis not in model_name:
          error_str = "axis(%s) not in dm" % axis
          raise Exception(error_str)

    # model grid
    x = [np.array(dm[s]) for s in model_name]
    xx = np.meshgrid(*x, indexing='ij')
    nn = xx[0].size

    # calculate cc values over all grid points 
    dm_grid = {}
    for i in range(model_num):
      par = model_name[i]
      dm_grid[par] = xx[i].flatten()

    zz, weight = self.cc_perturbed_seismogram(dm=dm_grid, plot=plot_seism)
    zz /= weight # weighted average of normalized zero-lag CC
    zz = zz.reshape(xx[0].shape)

    # get maximum cc value 
    imax = np.argmax(zz)
    idx_max = np.unravel_index(imax, xx[0].shape)
    zz_max = zz[idx_max]

    # plot out results
    with PdfPages(outfig) as pdf:
      for xy in axes:
        # plot 1D section
        if (len(xy) == 1) or (xy[0] == xy[1]):
          ix = model_name.index(xy[0])
          idx = [range(len(v)) for v in x]
          for i in range(model_num):
            if i != ix:
              idx[i] = [idx_max[i],]
          zz_1d = zz[np.ix_(*idx)].squeeze()

          # plot cross-sections through maximum point
          #fig = plt.figure(figsize=(8.5,11))
          fig = plt.figure()
          plt.plot(x[ix], zz_1d, 'ro')

          xlabel_str = "alpha (* d%s)" % model_name[ix]
          plt.xlabel(xlabel_str)
          plt.ylabel("weighted avg. CC")

          pdf.savefig()
          plt.close()

        # plot 2D cross section
        else:
          ix = model_name.index(xy[0])
          iy = model_name.index(xy[1])
          idx = [range(len(v)) for v in x]
          for i in range(model_num):
            if i != ix and i != iy:
              idx[i] = [idx_max[i],]
          idx = np.ix_(*idx)
          zz_2d = zz[idx].squeeze()
          xx_2d = xx[ix][idx].squeeze()
          yy_2d = xx[iy][idx].squeeze()
          if ix > iy:
            zz_2d = np.transpose(zz_2d)
            xx_2d = np.transpose(xx_2d)
            yy_2d = np.transpose(yy_2d)
          # plot 2D surface through the maximum
          #fig = plt.figure(figsize=(8.5, 11))
          fig = plt.figure()

          title_str = "average weighted normalized zero-lag CC"
          plt.title(title_str)

          levels = zz_max * np.linspace(0.5, 1.0, 20)
          cs = plt.contour(xx_2d, yy_2d, zz_2d, levels, colors='k')
          plt.clabel(cs, fontsize=9, inline=1)

          x_max = xx[ix][idx_max]
          y_max = xx[iy][idx_max]
          text_str = "(%.1f,%.1f)" % (x_max, y_max)
          plt.text(x_max, y_max, text_str,
              horizontalalignment="center", verticalalignment="top")
          plt.plot(x_max, y_max, 'ro')

          xlabel_str = "alpha (* d%s)" % model_name[ix]
          plt.xlabel(xlabel_str)
          ylabel_str = "alpha (* d%s)" % model_name[iy]
          plt.ylabel(ylabel_str)

          pdf.savefig()
          plt.close()

    # get best model
    event = self.data['event']
    src_perturb = self.data['src_perturb']
    
    model = {}
    for par in dm:
      m0 = event[par]
      dm = src_perturb[par]
      ix = model_name.index(par)
      alpha = xx[ix][idx_max]
      model[par] = {
          'm0':m0, 
          'dm':dm,
          'alpha':alpha,
          'm1':m0+alpha*dm
          }

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
    if 'mt' in model: mt = model['mt']['m1']

    with open(cmt_file, 'w') as fp:
      fp.write('%s | wCCmax %f weight %f\n'
          % (' '.join(header), zz_max, weight))
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

    return model, zz_max, weight

#
#======================================================
#

  def relocate_1d(self, 
      event_id,
      window_id_list=['F.p,P', 'F.s,S'], 
      fix_depth=False,
      out_cmt_file=None):
    """relocate event using ray path in reference earth model
    """
    # check inputs
    events = self.data['events']
    if event_id not in events:
      raise Exception("[ERROR] %s does NOT exist. Exit" % (event_id))
      sys.exit()

    # select windows
    sta_win_id_list = []
    event = events[event_id]
    stations = event['stations']
    for station_id in stations:

      station = stations[station_id]
      if station['stat']['code'] < 0:
        continue

      windows = station['windows']
      for window_id in window_id_list:
        if window_id not in windows:
          continue

        window = windows[window_id]
        if window['stat']['code'] != 1:
          continue 

        misfit = window['misfit']
        #if window['quality']['SNR'] < min_SNR or \
        #    misfit['CC0'] < min_CC0 or \
        #    misfit['CCmax'] < min_CCmax or\
        #    abs(misfit['CC_time_shift']) > max_CC_time_shift:
        #  continue

        sta_win_id = (station_id, window_id)
        sta_win_id_list.append(sta_win_id)

    # create sensitivity matrix G in local NED coordinate
    # G * dm  = dt_cc
    # G: [[-px_1, -py_1, -pz_1, 1.0], # ray1
    #   [-px_2, -py_2, -pz_2, 1.0], # ray2
    #   ...]
    # dm: [dNorth(km), dEast, dDepth, dT(sec)]
    # dt_cc: [dt1, dt2, ...]
    n = len(sta_win_id_list)
    G = np.zeros((n, 4))
    dt_cc = np.zeros(n)
    R_Earth_km = 6371.0
    gcmt = event['gcmt']
    evdp = gcmt['depth']
    for i in range(n):
      sta_win_id = sta_win_id_list[i]
      station_id = sta_win_id[0]
      window_id = sta_win_id[1]
  
      station = stations[station_id]
      meta = station['meta']
      window = station['windows'][window_id]
      phase = window['phase']
      misfit = window['misfit']
      weight = window['weight']

      azimuth = np.deg2rad(meta['azimuth'])
      takeoff_angle = phase['takeoff_angle']
      takeoff = np.deg2rad(takeoff_angle + 180.0*(takeoff_angle<0))
      ray_param = phase['ray_param']
      slowness = ray_param / (R_Earth_km - evdp) #unit: s/km
      # local coordinate: NED
      pd = np.cos(takeoff) * slowness
      pn = np.cos(azimuth) * np.sin(takeoff) * slowness
      pe = np.sin(azimuth) * np.sin(takeoff) * slowness
      # create sensitivity matrix 
      G[i,:] = weight * np.array([-pn, -pe, -pd, 1.0]) # -p: from receiver to source
      dt_cc[i] = weight * misfit['CC_time_shift']

    #linearized inversion (can be extended to second order using dynamic ray-tracing)
    if fix_depth: 
      G[:, 2] = 0.0
    dm, residual, rank, sigval = np.linalg.lstsq(G, dt_cc)

    # convert dm from NED to ECEF coordinate
    evla = gcmt['latitude']
    evlo = gcmt['longitude']

    slat = np.sin(np.deg2rad(evla))
    clat = np.cos(np.deg2rad(evla))
    slon = np.sin(np.deg2rad(evlo))
    clon = np.cos(np.deg2rad(evlo))

    N = [-slat*clon, -slat*slon, clat]
    E = [-slon, clon, 0.0]
    D = [-clat*clon, -clat*slon, -slat]

    ev_dx = N[0]*dm[0] + E[0]*dm[1] + D[0]*dm[2]
    ev_dy = N[1]*dm[0] + E[1]*dm[1] + D[1]*dm[2]
    ev_dz = N[2]*dm[0] + E[2]*dm[1] + D[2]*dm[2]
    ev_dt = dm[3]

    # initialize pyproj objects
    geod = pyproj.Geod(ellps='WGS84')
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    # old location in ECEF (meters)
    evx, evy, evz = pyproj.transform(lla, ecef, evlo, evla, -1000.0*evdp)

    # new location in ECEF (meters)
    evx1 = evx + ev_dx*1000.0
    evy1 = evy + ev_dy*1000.0
    evz1 = evz + ev_dz*1000.0
    # in LLA
    evlo1, evla1, evalt1 = pyproj.transform(ecef, lla, evx1, evy1, evz1)
    evdp1 = -evalt1/1000.0

    # residuals 
    # linearized modelling
    dt_syn = G.dot(dm)
    dt_res = dt_cc - dt_syn

    # make results
    new_centroid_time = UTCDateTime(gcmt['centroid_time']) + ev_dt
    reloc_dict = {
        'window_id_list': window_id_list,
        'singular_value': sigval.tolist(),
        'dm': {'dNorth':dm[0], 'dEast':dm[1], 'dDepth':dm[2], 
          'dT':dm[3]},
        'latitude':evla1, 
        'longitude':evlo1, 
        'depth':evdp1,
        'centroid_time': str(new_centroid_time),
        'data': {'num':n, 'mean':np.mean(dt_cc), 'std':np.std(dt_cc)},
        'residual': {'mean':np.mean(dt_res), 'std':np.std(dt_res)} }

    event['relocate'] = reloc_dict

    # make new CMTSOLUTION file
    if out_cmt_file:
      M = gcmt['moment_tensor']
      with open(out_cmt_file, 'w') as fp:
        # header line: 
        #PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
        # which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region
        fp.write(new_centroid_time.strftime(
          'RELOC %Y %m %d %H %M %S.%f ') + \
          '%.4f %.4f %.1f 0.0 0.0 END\n' % (evla1,evlo1,evdp1) )
        fp.write('event name:    %s\n'   % (event_id))
        fp.write('time shift:    0.0\n'        ) 
        fp.write('tau:   %.1f\n'   % (gcmt['tau']))
        #fp.write('half duration:   0.0\n'  % (gcmt['tau']))
        fp.write('latitude:    %.4f\n'   % (evla1)   )
        fp.write('longitude:     %.4f\n'   % (evlo1)   )
        fp.write('depth:       %.4f\n'   % (evdp1)   )
        fp.write('Mrr:       %12.4e\n' % (M[0][0]) )
        fp.write('Mtt:       %12.4e\n' % (M[1][1]) )
        fp.write('Mpp:       %12.4e\n' % (M[2][2]) )
        fp.write('Mrt:       %12.4e\n' % (M[0][1]) )
        fp.write('Mrp:       %12.4e\n' % (M[0][2]) )
        fp.write('Mtp:       %12.4e\n' % (M[1][2]) )

#
#======================================================
#


  def plot_misfit(self, event_id, window_id, out_file=None):
    """Plot misfit for a certain event and window_id  
    """
    # CC0 map  | CC0 v.s. SNR (size ~ weight)
    #------------|-----------------
    # DTcc map   | avg. CC0      

    # check inputs
    events = self.data['events']
    if event_id not in events:
      raise Exception("[ERROR] %s does NOT exist. Exit" % (event_id))
      sys.exit()
    event = events[event_id]
    stations = event['stations']

    # get list of station,window id
    #sta_win_id_list = []
    stla_list = []
    stlo_list = []
    cc_dt_list = []
    CC0_list = []
    CCmax_list = []
    snr_list = []
    weight_list = []
    for station_id in stations:
      station = stations[station_id]
      windows = station['windows']

      # skip bad station 
      if station['stat']['code'] < 0:
        continue

      if window_id not in windows:
        continue

      window = windows[window_id]
      if window['stat']['code'] != 1:
        continue

      meta = station['meta']
      misfit = window['misfit']
      quality = window['quality']

      #sta_win_id = (station_id, window_id)
      #sta_win_id_list.append(sta_win_id)
      stla_list.append(meta['latitude'])
      stlo_list.append(meta['longitude'])
      cc_dt_list.append(misfit['CC_time_shift'])
      CC0_list.append(misfit['CC0'])
      CCmax_list.append(misfit['CCmax'])
      snr_list.append(quality['SNR'])
      weight_list.append(window['weight'])

    # get event data
    gcmt = event['gcmt']
    evla = gcmt['latitude']
    evlo = gcmt['longitude']
    M = gcmt['moment_tensor']
    Mrr = M[0][0]
    Mtt = M[1][1]
    Mpp = M[2][2]
    Mrt = M[0][1]
    Mrp = M[0][2]
    Mtp = M[1][2]
    focmec = [ Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]

    # map range
    min_lat = min(min(stla_list), evla)
    max_lat = max(max(stla_list), evla)
    lat_range = max_lat - min_lat
    min_lat -= 0.1*lat_range
    max_lat += 0.1*lat_range
    min_lon = min(min(stlo_list), evlo)
    max_lon = max(max(stlo_list), evlo)
    lon_range = max_lon - min_lon
    min_lon -= 0.1*lon_range
    max_lon += 0.1*lon_range
    #lat_true_scale = np.mean(stla_list)
    lat_0 = np.mean(stla_list)
    lon_0 = np.mean(stlo_list)
    # 
    parallels = np.arange(0.,81,10.)
    meridians = np.arange(0.,351,10.)

    # figure size
    fig = plt.figure(figsize=(11, 8.5))
    str_title = '%s %s' % (event_id, window_id)
    fig.text(0.5, 0.95, str_title, size='x-large', 
        horizontalalignment='center')

    #------ color map CC_time_shift, symbol size ~ SNR 
    ax = fig.add_axes([0.05, 0.5, 0.4, 0.35])
    ax.set_title("DT_cc (symbol_size ~ SNR)")

    m = Basemap(projection='merc', resolution='l',
        llcrnrlat=min_lat, llcrnrlon=min_lon, 
        urcrnrlat=max_lat, urcrnrlon=max_lon,
        lat_0=lat_0, lon_0=lon_0 )
    m.drawcoastlines(linewidth=0.1)
    m.drawcountries(linewidth=0.1)
    m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1])
    m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1])
    
    # CC_time_shift, SNR
    sx, sy = m(stlo_list, stla_list)
    size_list = [ 0.1 if x<0.1 else x for x in snr_list ]
    im = m.scatter(sx, sy, s=size_list, marker='o',
        c=cc_dt_list, cmap='seismic', 
        edgecolor='grey', linewidth=0.05)
    mean_amp = np.mean(cc_dt_list)
    std_amp = np.std(cc_dt_list)
    #plot_amp = abs(mean_amp)+std_amp
    plot_amp = 5.0 
    im.set_clim(-plot_amp, plot_amp)
    
    # focal mechanism
    sx, sy = m(evlo, evla)
    b = Beach(focmec, xy=(sx, sy), width=200000, linewidth=0.2, 
        facecolor='k')
    ax.add_collection(b)
    
    # colorbar
    cbar_ax = fig.add_axes([0.46, 0.575, 0.005, 0.2])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cbar_ax.tick_params(labelsize=9)
    cbar_ax.set_xlabel('DT_cc(s)', fontsize=9)
    cbar_ax.xaxis.set_label_position('top')
   
    #------ color map CC0, symbol size ~ SNR
    ax = fig.add_axes([0.05, 0.05, 0.4, 0.35])
    ax.set_title("CC0 (symbol_size ~ SNR)")

    m = Basemap(projection='merc', resolution='l',
        llcrnrlat=min_lat, llcrnrlon=min_lon, 
        urcrnrlat=max_lat, urcrnrlon=max_lon,
        lat_0=lat_0, lon_0=lon_0 )
    m.drawcoastlines(linewidth=0.1)
    m.drawcountries(linewidth=0.1)
    m.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1])
    m.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1])
    
    # CC0, SNR 
    sx, sy = m(stlo_list, stla_list)
    #size_list = [ 20**x for x in CCmax_list ]
    size_list = [ 0.1 if x<0.1 else x for x in snr_list ]
    im = m.scatter(sx, sy, s=size_list, marker='o',
        c=CC0_list, cmap='jet', 
        edgecolor='grey', linewidth=0.05)
    im.set_clim(0.5, 1.0)
    
    # focal mechanism
    sx, sy = m(evlo, evla)
    b = Beach(focmec, xy=(sx, sy), width=200000, linewidth=0.2, 
        facecolor='k')
    ax.add_collection(b)
 
    #add colorbar
    cbar_ax = fig.add_axes([0.46, 0.125, 0.005, 0.2])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cbar_ax.tick_params(labelsize=9) 
    cbar_ax.set_xlabel('CC0', fontsize=9)
    cbar_ax.xaxis.set_label_position('top')

    #------ SNR v.s. CC0, colored by cc_dt, size ~ weight
    ax = fig.add_axes([0.58, 0.65, 0.35, 0.2])
    im = ax.scatter(snr_list, CC0_list, marker='o', 
        s=10.*np.array(weight_list), 
        c=cc_dt_list, cmap='seismic',
        edgecolor='grey', linewidth=0.05)
    mean_amp = np.mean(cc_dt_list)
    std_amp = np.std(cc_dt_list)
    #plot_amp = abs(mean_amp)+std_amp
    plot_amp = 5.0
    im.set_clim(-plot_amp, plot_amp)
    #ax.set_xlim([min(snr_list), max(snr_list)])
    #ax.set_ylim([min(CCmax_list), max(CCmax_list)])
    ax.set_xlim([0, max(snr_list)])
    ax.set_ylim([0.3, 1.0])
    ax.set_xlabel("SNR")
    ax.set_ylabel("CC0")
    #add colorbar
    cbar_ax = fig.add_axes([0.95, 0.65, 0.005, 0.2])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cbar_ax.tick_params(labelsize=9)
    cbar_ax.set_xlabel('DT_cc(s)', fontsize=9)
    cbar_ax.xaxis.set_label_position('top')

    ##------ CC0 v.s. CCmax, colored by cc_dt
    #ax = fig.add_axes([0.58, 0.375, 0.35, 0.2])
    #im = ax.scatter(CC0_list, CCmax_list, marker='o', s=10,
    #    c=cc_dt_list, cmap='seismic',
    #    edgecolor='grey', linewidth=0.05)
    #mean_amp = np.mean(cc_dt_list)
    #std_amp = np.std(cc_dt_list)
    #plot_amp = abs(mean_amp)+std_amp
    #im.set_clim(-plot_amp, plot_amp)
    #ax.set_xlim([min(CC0_list), max(CC0_list)])
    #ax.set_ylim([min(CCmax_list), max(CCmax_list)])
    #ax.set_xlabel("CC0")
    #ax.set_ylabel("CCmax")
    ##add colorbar
    #cbar_ax = fig.add_axes([0.95, 0.375, 0.005, 0.2])
    #fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    #cbar_ax.tick_params(labelsize=9)
    #cbar_ax.set_xlabel('cc_dt(s)', fontsize=9)
    #cbar_ax.xaxis.set_label_position('top')

    ##------ cc_dt v.s. CCmax, colored by SNR
    #ax = fig.add_axes([0.58, 0.1, 0.35, 0.2])
    #im = ax.scatter(cc_dt_list, CCmax_list, marker='o', s=10, 
    #    c=snr_list, cmap='seismic',
    #    edgecolor='grey', linewidth=0.05)
    #im.set_clim(min(snr_list), max(snr_list))
    #ax.set_xlim([min(cc_dt_list), max(cc_dt_list)])
    #ax.set_ylim([min(CCmax_list), max(CCmax_list)])
    #ax.set_xlabel("cc_dt")
    #ax.set_ylabel("CCmax")
    ##add colorbar
    #cbar_ax = fig.add_axes([0.95, 0.1, 0.005, 0.2])
    #fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    #cbar_ax.tick_params(labelsize=9)
    #cbar_ax.set_xlabel('SNR(dB)', fontsize=9)
    #cbar_ax.xaxis.set_label_position('top')

    ##------ histogram of dt_cc and dt_res
    #ax1 = fig.add_axes([0.5,0.28,0.4,0.15])
    #n, bins, patches = ax1.hist(dt_cc, 50, facecolor='green', alpha=0.75)
    #amp = max(abs(dt_cc))
    #ax1.set_xlim([-amp, amp])
    #ax1.set_title('dt_cc: mean=%.2f std=%.2f' % (np.mean(dt_cc), np.std(dt_cc)))
    #ax1.tick_params(labelsize=10) 
    #
    #ax2 = fig.add_axes([0.5,0.07,0.4,0.15])
    #n, bins, patches = ax2.hist(dt_res, 50, facecolor='green', alpha=0.75)
    #amp = max(abs(dt_cc))
    #ax2.set_xlim([-amp, amp])
    #ax2.set_title('dt_res: mean=%.2f std=%.2f' % (np.mean(dt_res), np.std(dt_res)))
    #ax2.set_xlabel('dt (sec)')
    #ax2.tick_params(labelsize=10)

    #------ save figure
    if not out_file:
      out_file = '%s_%s.pdf' % (event_id, window_id)
    fig.savefig(out_file, format='pdf')
    #fig.savefig("misfit.pdf", bbox_inches='tight', format='pdf')

#
#======================================================
#

  def plot_seismogram(self,
      savefig=False, out_dir='plot',
      plot_param={
        'time':[0,100], 'rayp':10., 'azbin':10, 'window_id':'F.p,P',
        'SNR':None, 'CC0':None, 'CCmax':None, 'dist':None }
      ):
    """ Plot seismograms for one event
      azbin: azimuthal bin size
      win: 
    """
    comp_name = ['R', 'T', 'Z']
    sampling_rate = self.data['sampling_rate']
    data_delta = 1.0/sampling_rate
    data_nyq = 0.5*sampling_rate

    #------ selection parameters
    plot_time = plot_param['time']
    plot_azbin = plot_param['azbin']
    plot_rayp = plot_param['rayp']
    plot_window_id = plot_param['window_id']

    plot_SNR = plot_param['SNR']
    plot_CC0 = plot_param['CC0']
    plot_CCmax = plot_param['CCmax']
    plot_dist = plot_param['dist']

    #------ event info
    event = self.data['event']
    t0 = event['t0']
    tau = event['tau']
    evla = event['latitude']
    evlo = event['longitude']
    evdp = event['depth']
    mt = event['mt_rtp']
    Mrr = mt[0][0]
    Mtt = mt[1][1]
    Mpp = mt[2][2]
    Mrt = mt[0][1]
    Mrp = mt[0][2]
    Mtp = mt[1][2]
    focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]

    #------ station info
    station_dict = self.data['station']
    stla_all = []
    stlo_all = []
    dist_all = []
    for station_id in station_dict:
      station = station_dict[station_id]
      meta = station['meta']
      window_dict = station['window']
      # select data 
      if station['stat']['code'] < 0:
        continue
      if plot_window_id not in window_dict:
        continue
      stla_all.append(meta['latitude'])
      stlo_all.append(meta['longitude'])
      dist_all.append(meta['dist_degree'])

    #------ traveltime curve
    model = TauPyModel(model="ak135")
    dist_ttcurve = np.arange(0.0,max(dist_all),0.5)
    ttcurve_p = []
    ttcurve_P = []
    ttcurve_s = []
    ttcurve_S = []
    for dist in dist_ttcurve:
      arrivals = model.get_travel_times(
          source_depth_in_km=evdp, 
          distance_in_degree=dist, 
          phase_list=['p','P','s','S'])
      for arr in arrivals:
        if arr.name == 'p':
          ttcurve_p.append((arr.distance, arr.time, arr.ray_param))
        elif arr.name == 'P':
          ttcurve_P.append((arr.distance, arr.time, arr.ray_param))
        elif arr.name == 's':
          ttcurve_s.append((arr.distance, arr.time, arr.ray_param))
        elif arr.name == 'S':
          ttcurve_S.append((arr.distance, arr.time, arr.ray_param))
    # sort phases
    ttcurve_p = sorted(ttcurve_p, key=lambda x: x[2])
    ttcurve_P = sorted(ttcurve_P, key=lambda x: x[2])
    ttcurve_s = sorted(ttcurve_s, key=lambda x: x[2])
    ttcurve_S = sorted(ttcurve_S, key=lambda x: x[2])

    #------ map configuration 
    min_lat = min(min(stla_all), evla)
    max_lat = max(max(stla_all), evla)
    lat_range = max_lat - min_lat
    min_lat -= 0.1*lat_range
    max_lat += 0.1*lat_range
    min_lon = min(min(stlo_all), evlo)
    max_lon = max(max(stlo_all), evlo)
    lon_range = max_lon - min_lon
    min_lon -= 0.1*lon_range
    max_lon += 0.1*lon_range
    lat_0 = np.mean(stla_all)
    lon_0 = np.mean(stlo_all)
    # 
    parallels = np.arange(0.,81,10.)
    meridians = np.arange(0.,351,10.)

    #------ plot azimuthal bins (one figure per azbin)
    if plot_azbin <= 0.0:
      raise Exception("plot_param['azbin']=%f must > 0.0" % plot_azbin)

    for az in np.arange(0, 360, plot_azbin):
      azmin = az
      azmax = az + plot_azbin

      #print(azmin, azmax

      #---- gather data for the current azbin 
      data_azbin = {}
      for station_id in station_dict:
        # skip bad station
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue

        # skip un-selected station 
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if plot_dist:
          if dist_degree < min(plot_dist) or dist_degree > max(plot_dist):
            continue
        if azimuth < azmin or azimuth >= azmax:
          continue
        if plot_window_id not in window_dict:
          continue

        # skip bad window
        window_dict = station['window']
        window = window_dict[plot_window_id]
        quality = window['quality']
        cc = window['cc']
        if window['stat']['code'] <= 0:
          continue
        if plot_SNR and quality['SNR']<min(plot_SNR):
          continue
        if plot_CC0 and cc['CC0']<min(plot_CC0):
          continue
        if plot_CCmax and cc['CCmax']<min(plot_CCmax):
          continue

        # get seismograms: syn/obs
        waveform = station['waveform']
        data_starttime = waveform['starttime']
        data_npts = waveform['npts']
        obs = waveform['obs']
        grf = waveform['grf']
        # filter parameter
        filter_param = window['filter']
        filter_a = filter_param['a']
        filter_b = filter_param['b']
        # filter seismograms 
        obs = signal.lfilter(filter_b, filter_a, obs)
        grf = signal.lfilter(filter_b, filter_a, grf)
        # convolve stf on grf
        freq = np.fft.rfftfreq(data_npts, d=data_delta)
        src_spectrum = stf_gauss_spectrum(freq, event['tau'])
        syn = np.fft.irfft(src_spectrum*np.fft.rfft(grf), data_npts)
        # project to polarity defined by the window
        proj_matrix = window['polarity']['proj_matrix']
        obs = np.dot(proj_matrix, obs)
        syn = np.dot(proj_matrix, syn)
        # rotate EN(0:2) -> RT(0:2) (T-R-Z: right-hand convention)
        Raz = (meta['back_azimuth'] + 180.0) % 360.0
        sin_Raz = np.sin(np.deg2rad(Raz))
        cos_Raz = np.cos(np.deg2rad(Raz))
        proj_matrix = [ 
            [sin_Raz,  cos_Raz], 
            [cos_Raz, -sin_Raz] ]
        obs[0:2,:] = np.dot(proj_matrix, obs[0:2,:])
        syn[0:2,:] = np.dot(proj_matrix, syn[0:2,:])
        # append to data
        data_dict = {
            'meta': meta,
            'window': window,
            'syn': syn,
            'obs': obs
            }
        data_azbin[station_id] = data_dict
      #endfor station_id in station_dict:
  
      #---- skip empty azbin
      if not data_azbin:
        warn_str = "No station in the azbin [%f %f]." %(azmin, azmax)
        warnings.warn(warn_str)
        continue

      #---- create figure
      fig = plt.figure(figsize=(8.5, 11)) # US letter
      str_title = '{:s} (win: {:s}, az: {:04.1f}~{:04.1f})'.format(
          event['id'], plot_window_id, azmin, azmax)
      fig.text(0.5, 0.95, str_title, size='x-large', horizontalalignment='center')
      #---- plot station/event map
      ax_origin = [0.3, 0.74]
      ax_size = [0.4, 0.2]
      ax_map = fig.add_axes(ax_origin + ax_size)
      ax_bm = Basemap(projection='merc', resolution='l',
          llcrnrlat=min_lat, llcrnrlon=min_lon, 
          urcrnrlat=max_lat, urcrnrlon=max_lon,
          lat_0=lat_0, lon_0=lon_0 )
      ax_bm.drawcoastlines(linewidth=0.1)
      ax_bm.drawcountries(linewidth=0.1)
      ax_bm.drawparallels(parallels, linewidth=0.1, labels=[1,0,0,1], 
          fontsize=10, fmt='%3.0f')
      ax_bm.drawmeridians(meridians, linewidth=0.1, labels=[1,0,0,1], 
          fontsize=10, fmt='%3.0f')
      sx, sy = ax_bm(stlo_all, stla_all)
      ax_bm.scatter(sx, sy, s=10, marker='^', facecolor='blue', edgecolor='')
      # plot focal mechanism
      sx, sy = ax_bm(evlo, evla)
      bb_width = 110000.0 * np.abs(max(stlo_all)-min(stlo_all)) * 0.1
      b = Beach(focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor='r')
      ax_map.add_collection(b)
      #-- plot the station location
      stla = [ x['meta']['latitude'] for x in data_azbin.itervalues() ]
      stlo = [ x['meta']['longitude'] for x in data_azbin.itervalues() ]
      sx, sy = ax_bm(stlo, stla)
      ax_bm.scatter(sx, sy, s=10, marker='^', facecolor='red', edgecolor='')

      #-- create axis for seismograms
      ax_RTZ = []
      for i in range(3):
        ax_origin = [0.07+0.3*i, 0.05]
        ax_size = [0.25, 0.65]
        ax_RTZ.append(fig.add_axes(ax_origin + ax_size))

      #-- plot traveltime curves
      for i in range(3):
        ax = ax_RTZ[i]
        ax.plot([x[1]-plot_rayp*x[0] for x in ttcurve_p], \
            [x[0] for x in ttcurve_p], 'b-', linewidth=0.2)
        ax.plot([x[1]-plot_rayp*x[0] for x in ttcurve_P], \
            [x[0] for x in ttcurve_P], 'b-', linewidth=0.2)
        ax.plot([x[1]-plot_rayp*x[0] for x in ttcurve_s], \
            [x[0] for x in ttcurve_s], 'c-', linewidth=0.2)
        ax.plot([x[1]-plot_rayp*x[0] for x in ttcurve_S], \
            [x[0] for x in ttcurve_S], 'c-', linewidth=0.2)

      #-- ylim setting
      y = [ x['meta']['dist_degree'] for x in data_azbin.itervalues() ]
      ny = len(y)
      plot_dy = 0.5*(max(y)-min(y)+1)/ny
      if plot_dist:
        plot_ymax = max(plot_dist) + 2*plot_dy
        plot_ymin = min(plot_dist) - 2*plot_dy
      else:
        plot_ymax = max(y) + 2*plot_dy
        plot_ymin = min(y) - 2*plot_dy

      #-- plot each station
      for station_id in data_azbin:
        station = data_azbin[station_id]
        meta = station['meta']
        window = station['window']
        syn = station['syn']
        obs = station['obs']

        # get plot time 
        dist_degree = meta['dist_degree']
        reduced_time = dist_degree * plot_rayp
        # time of first sample referred to centroid time 
        t0 = data_starttime - event['t0']
        # time of samples referred to centroid time
        syn_times = data_delta*np.arange(data_npts) + t0
        # plot time window
        plot_t0 = min(plot_time) + reduced_time
        plot_t1 = max(plot_time) + reduced_time
        plot_idx = (syn_times > plot_t0) & (syn_times < plot_t1)
        # plot time
        t_plot = syn_times[plot_idx] - reduced_time

        #  window begin/end
        taper = window['taper']
        win_starttime = taper['starttime']
        win_endtime = taper['endtime']
        win_t0 = (win_starttime - event['t0']) - reduced_time
        win_t1 = (win_endtime - event['t0']) - reduced_time

        # plot seismograms
        Amax_obs = np.sqrt(np.max(np.sum(obs[:,plot_idx]**2, axis=0)))
        Amax_syn = np.sqrt(np.max(np.sum(syn[:,plot_idx]**2, axis=0)))
        for i in range(3):
          ax = ax_RTZ[i]
          ax.plot(t_plot, plot_dy*obs[i,plot_idx]/Amax_obs+dist_degree, \
              'k-', linewidth=0.5)
          ax.plot(t_plot, plot_dy*syn[i,plot_idx]/Amax_syn+dist_degree, \
              'r-', linewidth=0.5)
          # mark measure window range
          ax.plot(win_t0, dist_degree, 'k|', markersize=8)
          ax.plot(win_t1, dist_degree, 'k|', markersize=8)
          # annotate amplitude
          if i == 0:
            ax.text(max(plot_time), dist_degree, '%.1e ' % (Amax_obs), 
                verticalalignment='bottom', 
                horizontalalignment='right', 
                fontsize=7, color='black')
            ax.text(max(plot_time), dist_degree, '%.1e ' % (Amax_syn), 
                verticalalignment='top', 
                horizontalalignment='right', 
                fontsize=7, color='red')
          # annotate CC0
          if i == 0:
            ax.text(max(plot_time), dist_degree, ' %.3f'%(window['cc']['CC0']),
                verticalalignment='center', fontsize=7)
          # annotate window weight
          if i == 1:
            ax.text(max(plot_time), dist_degree, ' %.1f' % (window['weight']),
                verticalalignment='center', fontsize=7)
          #annotate station names 
          if i == 2:
            #str_annot = '%.3f,%.1f,%s' % (
            #    misfit['CC0'], window['weight'], station_id)
            ax.text(max(plot_time), dist_degree, ' '+station_id, \
                verticalalignment='center', fontsize=7)
      #endfor data in data_azbin:

      #-- set axes limits and lables, annotation
      for i in range(3):
        ax = ax_RTZ[i]
        ax.set_xlim(min(plot_time), max(plot_time))
        ax.set_ylim(plot_ymin, plot_ymax)
        ax.set_title(comp_name[i])
        ax.set_xlabel('t - {:.1f}*dist (s)'.format(plot_rayp))
        ax.tick_params(axis='both',labelsize=10)
        # ylabel 
        if i == 0:
          ax.set_ylabel('dist (deg)')
        else:
          ax.set_yticklabels([])

      #-- save figures
      if savefig: 
        out_file = '%s/%s_az_%03d_%03d_%s.pdf' \
            % (out_dir, event['id'], azmin, azmax, plot_window_id)
        plt.savefig(out_file, format='pdf')
      else:
        plt.show()
      plt.close(fig)
    
#END class misfit