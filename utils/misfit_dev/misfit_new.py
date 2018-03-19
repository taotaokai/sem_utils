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
from event import Event
from station import Station 

#NOTE
# 1. spectrum relation between DFT and FT
#   x(n*dt): DFT[x]*dt ~ FT[x], IDFT[FT[x]]/dt ~ x

#====== utility functions
def is_equal(lst):
  return not lst or [lst[0]]*len(lst) == lst

def stf_gauss_spectrum(f, tau):
  """ spectrum of the Gaussian STF of unit area: 
      stf(t,tau) = 1/sqrt(PI)/tau * exp(-(t/tau)^2)
      F_stf = exp(- pi^2 * f^2 * tau^2)
  """
  F_src = np.exp(-np.pi**2 * f**2 * tau**2)
  return F_src

#======
class Misfit(object):
  """
  Class managing all misfit windows

  Unit: kg,m,s

  Coordinate: ECEF cartesian

  Attributes
  ----------
  event (class Event)

  stations[] (class Station) {
    status, history,
    network, station, location,
    metadata, 
    waveform,
    green, dgreen,
    windows[] (class Window),
  }

  src_frechet: {'t0':, 'tau':, 'xs':, 'mt':},
  src_perturb: {'t0':1.0, 'tau':1.0, 'xs':, 'mt':, },

  chi: # weighted average of normalized zero-lag CC

  }

  Methods
  -------
  setup_event
  setup_station
  read_observed_data:
      - seismograms are resampled to sampling_rate, with 0.8*nyq low-pass filter
  read_sythetic_green:
      - 0.8*nyq low-pass filter before resampling
  setup_window
  measure_window
  make_adjoint_source
  window_quality_control (determine bad/OK, and window weight)
  make_cmt_dxs/dmt
  waveform_der_dxs/dmt
  cc_perturbed_seisomgram
  grid_cc

  NOTE
  ----
  1D Earth model: ak135

  """
#======================================================
  def __init__(self):
    self.event = Event()
    self.stations = []
    self.misift_value = 0.0

#======================================================
  def save(self, out_file):
    """
    Use netCDF4 format for data persistency
    """
    from netCDF4 import Dataset
    root = Dataset(out_file, 'w')

    #------ event
    event = root.createGroup("/event")

    #------ station
    stations = root.createGroup("/station")

    for station in self.stations:
      station.


#======================================================
  def read_cmtsolution(self, cmt_file, isECEF=True):
    """
    Read in CMTSOLUTION
    """
    self.event.read_cmtsolution(cmt_file, isECEF=True)

#======================================================
  def add_stations(self, station_file,):
    """ 
    Add stations 

    Parameters
    ----------
    station_file : str 
      text file of net sta loc

    """
    #------ read station file
    with open(station_file, 'r') as f:
      lines = [x.replace('\n','').split('|')  \
          for x in f.readlines() if not(x.startswith('#'))]

    #------ add stations
    stid_list = [(x[0],x[1],x[2]) for x in lines]
    for sta in self.stations:
      stid = (sta.network, sta.station, sta.location)
      if stid in stid_list:
        sta.reset(net,sta,loc)
      else:
        self.stations.append(Station(net,sta,loc))

#======================================================
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

#======================================================
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

#======================================================
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

#======================================================
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

#======================================================
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

#======================================================
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