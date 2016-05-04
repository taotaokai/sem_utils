#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import warnings
#import os.path
import re
#
import numpy as np
import scipy.signal as signal
from scipy import interpolate
#
import pickle
#
from obspy import UTCDateTime, read, Trace
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
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


#====== utility functions
def is_equal(lst):
  return not lst or [lst[0]]*len(lst) == lst

def stf_gauss_spectrum(f, tau):
  """ spectrum of the Gaussian STF of unit area: 
      stf(t,tau) = 1/sqrt(PI)/tau * exp(-(t/tau)^2)
      F_stf = exp(- pi^2 * f^2 * tau^2)
  """
  F_src = np.exp(-np.pi**2 * f**2 * tau**2)
# F_ds_dt0 = -2.0j * np.pi * f * F_src
# F_ds_dtau = -2.0 * (np.pi*f)**2 * tau * F_src
# return F_src, F_ds_dt0, F_ds_dtau
  return F_src

#def stf_gauss(n, dt, tau):
#  """ Gaussian source time function
#    stf(t,tau) = 1/sqrt(PI)/tau * exp(-(t/tau)^2)
#    t = dt * [0:(n+1)/2, -n/2:-1]
#  """
#  n = int(n)
#  idx = np.arange(n)
#  idx[(n+1)/2:] -= n
#  t = idx * dt / tau
#  stf = np.exp(-t**2)/tau/np.pi**0.5
#  ds_dt0 = 2.0/tau * t * stf
#  ds_dtau = 2.0/tau * (t**2 - 0.5) * stf 
#  return stf, ds_dt0, ds_dtau

#======
class Misfit(object):
  """Class managing all misfit windows

  Unit: kg,m,s

  Coordinate: ECEF cartesian

  Data:
  {
    'event': {
        'stat':{'code':, 'msg':},
        'id':, # event ID
        'header':, # have centroid time in it
        'lattidue':, 'longitude':, 'depth':, 
        't0':[utcdatetime], # centroid time
        'tau':, # gaussian width exp(-(t-t0)^2/tau^2)/tau/pi^0.5
        'xs': [x, y, z], # source location in ECEF coordinate
        'mt':, # moment tensor in ECEF coord
        'mt_rtp':}, # moment tensor in spherical coord

    'src_frechet': { #derivative of misfit function w.r.t. source param.
        'stat':{'code':, 'msg':},
        't0':, 'dau':, 'xs':, 'mt':},

    'src_perturb': {'t0':1.0, 'tau':1.0, 'xs':, 'mt':, },

    'chi': # weighted average of normalized zero-lag CC

    'station': {

        <station_id> : {
            'stat': {code:, msg:},
            'meta': {
                latitude:, longitude:, elevation:, depth:,
                azimuth:, back_azimuth:, dist_degree:,
                'channel': [ {code,az,dip,...}, ...]
                'ttime':
            },

            'waveform': {
                'time_sample': {starttime:, delta:, nt:, nl:, nr},
                'obs': array([3,nt]), # observed seismograms
                'grf': array([3,nt], # green's function
                #'syn': 3 x nt (u), # synthetic seismograms (grf convolve stf)
            },

            'window': {
                <window_id>: {
                    'stat': {code:, msg:}, 
                    #bad:<0, ok:0, measured adj:1, and hessian_src:2

                    'filter': {'type':,'order':,'freqlim':,'a':,'b':},
                    'taper': {
                        'starttime':, 'endtime':,
                        'type':, 'ratio':, 'win':array(nt) }, 
                    'polarity': {
                        'component':, 'azimuth':, 'dip':,
                        'proj_matrix':array([3,3]) },

                    'quality': {
                        'Amax_obs':, 'Amax_noise':, 'Amax_syn':, 'SNR':, },
                    'cc': {'time':, 'cc':, 'cc_tshift': 
                        'CC0':, 'CCmax':, 'AR0':, 'ARmax':},
                    'weight':,

                    #'dchi_du':, # adjoint source from this window
                    #'dchi_dg':,

                    'hessian_src':{ # approximated Hessian of chi to source
                        ('dt0','dt0'):,
                        ('dtau','dtau'):,
                        ('dxs','dxs'):,
                        ...
                    },
                },
                <window_id>...
            },

            #'dchi_du': array([3,nt]), #chi: misfit function
            'dchi_dg': array([3,nt]), #conj(stf)*dchi_du

            'waveform_der': {
               'xs': {'dm':array(3), 'dg':, }, #finite-difference
               'mt': {'dm':array([3,3]), 'dg':,}, #linear in moment-tensor
            },

        },

        <station_id>...
    },
  }

  Methods:
    1. setup_event
    2. setup_station
    3. setup_window
    4. measure_window
    5. window_quality_control (determine bad/OK, and window weight)
    6. output_adj
    7. make_cmt_dxs/dmt
    8. waveform_der_dxs/dmt
    9. cc_perturbed_seisomgram
    10. grid_cc

  NOTE:
    1. 1D Earth model: ak135

  """
#
#======================================================
#

  def __init__(self):
    """ initialize
    """
    self.data = {
        'event':{},
        'src_frechet':{},
        'src_perturb':{'t0':1.0, 'tau':1.0},
        'station':{},
        'chi':0.0,
        }
#
#======================================================
#

  def save(self, filename='misfit.pkl'):
    """Save data
    """
    # use highest protocol available 
    with open(filename, 'wb') as fp:
      pickle.dump(self.data, fp, -1)
#
#======================================================
#

  def load(self, filename='misfit.pkl'):
    """Load data
    """
    with open(filename, 'rb') as fp:
      self.data = pickle.load(fp)
#
#======================================================
#

  def setup_event(self, cmt_file, ECEF=False):
    """cmt_file (str): CMTSOLUTION format file
    """
    with open(cmt_file, 'r') as f:
      lines = [ x for x in f.readlines() if not(x.startswith('#')) ]

    header = lines[0].split()
    year   = header[1]
    month  = header[2]
    day    = header[3]
    hour   = header[4]
    minute = header[5]
    second = header[6]

    lines = [x.split(":") for x in lines]
    event_id = lines[1][1].strip()
    time_shift = float(lines[2][1])

    # initialize pyproj objects
    geod = pyproj.Geod(ellps='WGS84')
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    if ECEF:
      tau = float(lines[3][1])
      x   = float(lines[4][1])
      y   = float(lines[5][1])
      z   = float(lines[6][1])
      # convert from ECEF(meters) to lla
      lon, lat, alt = pyproj.transform(ecef, lla, x, y, z)
      dep = -alt / 1000.0
    else:
      tau = float(lines[3][1]) / 1.628 # mimic triangle with gaussian
      lat = float(lines[4][1])
      lon = float(lines[5][1])
      dep = float(lines[6][1])
      # convert from lla to ECEF(meters)
      alt = -1000.0 * dep #NOTE ignore local topography
      x, y, z = pyproj.transform(lla, ecef, lon, lat, alt)

    # centroid time: t0
    isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(
        year, month, day, hour, minute, second)
    t0 = UTCDateTime(isotime) + time_shift
    # modify origin time in header line to have centroid time 
    header[1] = "{:04d}".format(t0.year)
    header[2] = "{:02d}".format(t0.month)
    header[3] = "{:02d}".format(t0.day)
    header[4] = "{:02d}".format(t0.hour)
    header[5] = "{:02d}".format(t0.minute)
    header[6] = "{:07.4f}".format(t0.second + 1.0e-6*t0.microsecond)

    # moment tensor
    # ECEF=false: 1,2,3 -> r,theta,phi
    # ECEF=true:  1,2,3 -> x,y,z
    m11 = float( lines[7][1])
    m22 = float( lines[8][1])
    m33 = float( lines[9][1])
    m12 = float(lines[10][1])
    m13 = float(lines[11][1])
    m23 = float(lines[12][1])
    mt = np.array([
      [m11, m12, m13], 
      [m12, m22, m23], 
      [m13, m23, m33]])
    # transform from spherical to cartesian coordinate
    r = (x**2 + y**2 + z**2)**0.5
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    # rotation matrix
    sthe = np.sin(theta)
    cthe = np.cos(theta)
    sphi = np.sin(phi)
    cphi = np.cos(phi)
    # basis transform matrix: e_x,y,z = a * e_r,t,p
    mt_xyz = np.zeros((3,3))
    mt_rtp = np.zeros((3,3))
    a = np.array(
        [ [ sthe*cphi, cthe*cphi, -1.0*sphi ],
          [ sthe*sphi, cthe*sphi,      cphi ],
          [ cthe     , -1.0*sthe,      0.0  ] ])
    if ECEF:
      mt_xyz = mt
      mt_rtp = np.dot(np.dot(np.transpose(a), mt), a)
    else: # spherical coordinate
      # harvard cmt use dyn*cm, change to N*m
      mt_rtp = mt*1.0e-7
      mt_xyz = np.dot(np.dot(a, mt_rtp), np.transpose(a))

    # add event
    event = {
        'id':event_id,
        'header':' '.join(header),
        'longitude':lon, 'latitude':lat, 'depth':dep,
        't0':t0, 'tau':tau,
        'xs':[x, y, z], 'mt':mt_xyz, 'mt_rtp':mt_rtp,
        'stat': {'code':0, 'msg':"created on "+UTCDateTime.now().isoformat()}
        }
    self.data['event'] = event

  #enddef setup_event

#
#======================================================
#

  def write_cmtsolution(self, out_file, ECEF=True):
    """ 
    Write out CMTSOLUTION file 
    """
    # modify origin time in header line to have centroid time 
    event = self.data['event']
    event_id = event['id']
    header = event['header']
    tau = event['tau']

    if ECEF:
      xs = event['xs']
      mt = event['mt']
      with open(out_file, 'w') as fp:
        fp.write('%s\n' % (header))
        fp.write('%-18s %s\n' % ('event name:', event_id))
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
    else:
      mt_rtp = event['mt_rtp']
      with open(out_file, 'w') as fp:
        fp.write('%s\n' % (header))
        fp.write('%-18s %s\n' % ('event name:', event_id))
        fp.write('%-18s %+15.8E\n' % ('time shift:',    0.0))
        fp.write('%-18s %+15.8E\n' % ('half duration:', tau*1.628))
        fp.write('%-18s %+15.8E\n' % ('latitude:',    event['latitude']))
        fp.write('%-18s %+15.8E\n' % ('longitude:',   event['longitude']))
        fp.write('%-18s %+15.8E\n' % ('depth:',       event['depth']))
        fp.write('%-18s %+15.8E\n' % ('Mrr(dyn*cm):', mt_rtp[0,0]*1.0e7))
        fp.write('%-18s %+15.8E\n' % ('Mtt(dyn*cm):', mt_rtp[1,1]*1.0e7))
        fp.write('%-18s %+15.8E\n' % ('Mpp(dyn*cm):', mt_rtp[2,2]*1.0e7))
        fp.write('%-18s %+15.8E\n' % ('Mrt(dyn*cm):', mt_rtp[0,1]*1.0e7))
        fp.write('%-18s %+15.8E\n' % ('Mrp(dyn*cm):', mt_rtp[0,2]*1.0e7))
        fp.write('%-18s %+15.8E\n' % ('Mtp(dyn*cm):', mt_rtp[1,2]*1.0e7))

#
#======================================================
#

  def setup_station(self,
      channel_file,
      band_code=None, 
      three_channels=True, 
      ):
    """ Setup station metadata.

      channel_file (str):
        FDSN-station text format file at channel level
      event_id_list (list): 
        list of event ID's to which stations are added [default: None]
        default to all events.
      band_code (str):
        instrument/band code [default: None]
      three_channels (bool):
        check for completeness of 3 channels [default: False]

      Note: 
      1) only use stations which have the same lat/lon/ele/depth 
        in all the available channels.
      2) gcmt info must be set first.
    """
    # initiate taup
    taup_model = TauPyModel(model="ak135")

    # read station file
    with open(channel_file, 'r') as f:
      lines = [x.replace('\n','').split('|')  \
          for x in f.readlines() if not(x.startswith('#'))]
    
    # get all station metadata
    metadata = {}
    for x in lines:
      net_sta_loc = (x[0], x[1], x[2])
      date1 = [ int(a) for a in re.sub("\D", " ", x[15]).split() ]
      date2 = [ int(a) for a in re.sub("\D", " ", x[16]).split() ]
      t1 = UTCDateTime(date1[0], date1[1], date1[2]) \
          + 60.0*(60.0*date1[3] + date1[4]) + date1[5]
      t2 = UTCDateTime(date2[0], date2[1], date2[2]) \
          + 60.0*(60.0*date2[3] + date2[4]) + date2[5]
      channel = {
          'code':    x[3],
          'latitude':  float(x[4]),
          'longitude':   float(x[5]),
          'elevation':   float(x[6]),
          'depth':     float(x[7]),
          'azimuth':   float(x[8]),
          'dip':     float(x[9]),
          'starttime':   t1,
          'endtime':   t2}
      if net_sta_loc not in metadata:
        metadata[net_sta_loc] = []
      metadata[net_sta_loc].append(channel)

    # check if event info is set 
    data = self.data
    if ('event' not in data) or (data['event']['stat']['code'] < 0):
      raise Exception('Event info not set.')
    event = data['event']

    # initialize station dict
    data['station'] = {}
    station = data['station']

    # station active time is set to centroid time
    # used to filter channel list
    active_time = event['t0']

    for net_sta_loc in metadata:
      # station_id: net.sta.loc
      station_id = '.'.join(net_sta_loc)

      # skip existing stations if not update
      if (station_id in station) and (not update):
        raise Exception('station %s already exist' % (station_id))

      # select channels which are active at the specified time 
      channel = [ x for x in metadata[net_sta_loc] 
          if x['starttime'] < active_time and x['endtime'] > active_time ]

      # select band code (e.g. BH )
      if band_code:
        n = len(band_code)
        channel = [ x for x in channel if x['code'][0:n]==band_code ]

      # check if all selected channels have the same location
      lats = [ x['latitude'] for x in channel ]
      lons = [ x['longitude'] for x in channel ]
      eles = [ x['elevation'] for x in channel ]
      deps = [ x['depth'] for x in channel ]
      if lats.count(lats[0])!=len(lats) or \
          lons.count(lons[0])!=len(lons) or \
          eles.count(eles[0])!=len(eles) or \
          deps.count(deps[0])!=len(deps):
        print("[WARNING] %s: " \
          "channels do NOT have the same coordinates, SKIP." \
          % (station_id))
        continue

      # check completeness for 3-components
      if three_channels:
        if len(channel) != 3:
          print('[WARNING] %s: not exactly 3 components found, '\
            'SKIP.' % (station_id))
          continue
        # check channel orientations
        Z_comp = [ (x['code'], x['azimuth'], x['dip'])
            for x in channel if x['code'][2] == 'Z']
        if len(Z_comp) != 1 or abs(Z_comp[0][2]) != 90.0: 
          print('[WARNING] %s: problematic Z channel, SKIP' \
              % (station_id))
          print('      channel: ', H_comp)
          continue
        H_comp = [ (x['code'], x['azimuth'], x['dip']) \
            for x in channel if x['code'][2] != 'Z']
        if len(H_comp) != 2 or \
            abs(H_comp[0][2]) != 0.0 or \
            abs(H_comp[1][2]) != 0.0 or \
            abs(np.cos(np.deg2rad(
              H_comp[0][1] - H_comp[1][1]))) > 0.1: 
          print('[WARNING] %s: problematic horizontal channels, SKIP'\
              % (station_id))
          print('      channel: ', H_comp)
          continue

      # geodetic and ak135 traveltimes
      dist, az, baz = gps2dist_azimuth(
          event['latitude'], event['longitude'],
          channel[0]['latitude'], channel[0]['longitude'])
      dist_degree = kilometer2degrees(dist/1000.0)

      arrivals = taup_model.get_travel_times(
          source_depth_in_km=event['depth'],
          distance_in_degree=dist_degree,
          phase_list=['ttp','tts'],
          )

      # make station metadata 
      meta = { #TODO remove dumplicated info in channels?
          'latitude': channel[0]['latitude'], 
          'longitude': channel[0]['longitude'],
          'elevation': channel[0]['elevation'],
          'depth': channel[0]['depth'],
          'channel': channel,
          'azimuth': az,
          'back_azimuth': baz,
          'dist_degree': dist_degree,
          'ttime': arrivals}

      # add station info
      station[station_id] = {
          'meta': meta,
          'waveform':{},
          'window':{},
          'waveform_der':{},
          'stat': {
            'code': 0,
            'msg': "created on "+UTCDateTime.now().isoformat()} 
          }

    #endfor net_sta_loc in stations_all:
  #enddef setup_stations_from_channel_file

#
#======================================================
#
  def radiation_pattern(self, 
      outfig="radiation.pdf", 
      phase_list=['P','S','pP','sS'],
      dist=[25, 180],
      ):
    """
    Plot radiation pattern in the local ENU coordinate

    Note:
    """
    #------ event location and focal mechanism 
    event = self.data['event']
    evla = event['latitude']
    evlo = event['longitude']
    evdp = event['depth']
    mt_rtp = event['mt_rtp']

    # locat ENU coord
    mt_enu = np.zeros((3,3))
    # e = p, n = -t, u = r
    mt_enu[0,0] = mt_rtp[2,2] # ee = pp
    mt_enu[0,1] = -1.0*mt_rtp[1,2] # en = -tp
    mt_enu[1,0] = mt_enu[0,1]
    mt_enu[0,2] = mt_rtp[0,2] # eu = rp
    mt_enu[2,0] = mt_enu[0,2]
    mt_enu[1,1] = mt_rtp[1,1] # nn = tt
    mt_enu[1,2] = mt_rtp[0,1] # nu = -rt
    mt_enu[2,1] = mt_enu[1,2]
    mt_enu[2,2] = mt_rtp[0,0] # uu = rr

    #------ stations 
    # initialize taup
    taup_model = TauPyModel(model="ak135")

    station_dict = self.data['station']
    for station_id in station_dict:
      station = station_dict[station_id]
      if station['stat']['code'] < 0:
        continue
      meta = station['meta']
      dist_degree = meta['dist_degree']
      if dist_degree > min(dist) and dist_degree < max(dist):
        arrivals = taup_model.get_travel_times(
            source_depth_in_km=evdp,
            distance_in_degree=dist_degree,
            phase_list=phase_list,
            )
        meta['ttime'] = arrivals

    #------ calculate radiation pattern for up-going rays 
    #-- outgoing rays 
    # cmpinc: incline angle measured from Up/Down direction for up/down-going rays
    # cmpaz: azimuth measured on EN plane clockwise from North
    dtheta = 0.5
    cmpinc = np.arange(0.0, 90.0+dtheta/2, dtheta) 
    cmpaz = np.arange(0.0, 360.0+dtheta/2, dtheta) 
    tt, pp = np.meshgrid(cmpinc, cmpaz)
    nn = tt.size
    sin_inc = np.sin(np.deg2rad(tt.flatten())) 
    cos_inc = np.cos(np.deg2rad(tt.flatten())) 
    sin_az = np.sin(np.deg2rad(pp.flatten())) 
    cos_az = np.cos(np.deg2rad(pp.flatten())) 
    #-- unit vectors of up-going rays 
    vu = np.zeros((3, nn))
    vu[0,:] = sin_inc * sin_az
    vu[1,:] = sin_inc * cos_az
    vu[2,:] = cos_inc
    #-- unit vectors of down-going rays 
    # inc is measured from Down
    vd = np.zeros((3, nn))
    vd[0,:] = sin_inc * sin_az
    vd[1,:] = sin_inc * cos_az
    vd[2,:] = -1.0*cos_inc
    #-- tractions of up/down-going rays
    mt_vu = np.dot(mt_enu, vu)
    mt_vd = np.dot(mt_enu, vd)
    #-- P radiation pattern
    Apu = np.zeros(nn) # amplitude of up-going rays
    Apd = np.zeros(nn) # amplitudes of down-going rays
    for i in range(nn):
      Apu[i] = np.dot(vu[:,i], mt_vu[:,i])
      Apd[i] = np.dot(vd[:,i], mt_vd[:,i])
    #-- S radiation pattern
    Ashu = np.zeros(nn)
    Ashd = np.zeros(nn)
    Asvu = np.zeros(nn)
    Asvd = np.zeros(nn)
    for i in range(nn):
      Fsu = mt_vu[:,i] - vu[:,i]*Apu[i]
      Fsd = mt_vd[:,i] - vd[:,i]*Apd[i]

      Ashu[i] = cos_az[i]*Fsu[0] - sin_az[i]*Fsu[1]
      Ashd[i] = cos_az[i]*Fsd[0] - sin_az[i]*Fsd[1]

      Asvu[i] = -1.0*cos_inc[i]*(sin_az[i]*Fsu[0] + cos_az[i]*Fsu[1]) \
          + sin_inc[i]*Fsu[2]
      Asvd[i] = 1.0*cos_inc[i]*(sin_az[i]*Fsd[0] + cos_az[i]*Fsd[1]) \
          + sin_inc[i]*Fsd[2]

    #-- reshape into 2D array
    Apu = Apu.reshape(tt.shape)
    Apd = Apd.reshape(tt.shape)
    Ashu = Ashu.reshape(tt.shape)
    Ashd = Ashd.reshape(tt.shape)
    Asvu = Asvu.reshape(tt.shape)
    Asvd = Asvd.reshape(tt.shape)

    #------ plot radiation pattern
    # polar to cartesian
    xx = tt * sin_az.reshape(tt.shape)
    yy = tt * cos_az.reshape(tt.shape)
    # outer circle
    az = np.arange(0.0, 360.1, 1)
    x_peri = 90.0 * np.sin(np.deg2rad(az))
    y_peri = 90.0 * np.cos(np.deg2rad(az))
    #
    with PdfPages(outfig) as pdf:

      #-- plot  P radiation pattern
      Ap_max = max(np.max(np.abs(Apu)), np.max(np.abs(Apd)))

      # up-going P
      fig = plt.figure()
      ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect='equal')
      # contour
      levels = np.linspace(-1.0,1.0,21)
      cs = ax.contour(xx, yy, Apu/Ap_max, levels, colors='k')
      ax.clabel(cs, fontsize='x-small', inline=1, fmt='%.1f')
      # plot axis
      ax.plot(x_peri, y_peri, 'k')
      for theta in np.arange(10, 90, 10):
        xc = theta * np.sin(np.deg2rad(az))
        yc = theta * np.cos(np.deg2rad(az))
        ax.plot(xc, yc, '-', color=(0.5,0.5,0.5), linewidth=0.5)
      for theta in np.arange(0, 360, 10):
        xc = 90.0 * np.sin(np.deg2rad(theta))
        yc = 90.0 * np.cos(np.deg2rad(theta))
        ax.plot([0, xc], [0, yc], '-', color=(0.5,0.5,0.5), linewidth=0.5)
      # plot stations
      for station_id in station_dict:
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if dist_degree > min(dist) and dist_degree < max(dist):
          for arr in meta['ttime']:
            if arr.name[0] == "p":
              takeoff_angle = arr.takeoff_angle
              inc = 180.0 - takeoff_angle
              xc = inc * np.sin(np.deg2rad(azimuth))
              yc = inc * np.cos(np.deg2rad(azimuth))
              ax.plot(xc, yc, 'r+')
      # set axis
      ax.set_title("Radiation pattern: up-going P")
      pdf.savefig()
      plt.close()

      # down-going P
      fig = plt.figure()
      ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect='equal')
      # contour
      levels = np.linspace(-1.0,1.0,21)
      cs = ax.contour(xx, yy, Apd/Ap_max, levels, colors='k')
      ax.clabel(cs, fontsize='x-small', inline=1, fmt='%.1f')
      # plot axis
      ax.plot(x_peri, y_peri, 'k')
      for theta in np.arange(10, 90, 10):
        xc = theta * np.sin(np.deg2rad(az))
        yc = theta * np.cos(np.deg2rad(az))
        ax.plot(xc, yc, '-', color=(0.5,0.5,0.5), linewidth=0.5)
      for theta in np.arange(0, 360, 10):
        xc = 90.0 * np.sin(np.deg2rad(theta))
        yc = 90.0 * np.cos(np.deg2rad(theta))
        ax.plot([0, xc], [0, yc], '-', color=(0.5,0.5,0.5), linewidth=0.5)
      # plot stations
      for station_id in station_dict:
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if dist_degree > min(dist) and dist_degree < max(dist):
          for arr in meta['ttime']:
            if arr.name[0] == "P":
              takeoff_angle = arr.takeoff_angle
              inc = takeoff_angle
              xc = inc * np.sin(np.deg2rad(azimuth))
              yc = inc * np.cos(np.deg2rad(azimuth))
              ax.plot(xc, yc, 'r+')
      # set axis
      ax.set_title("Radiation pattern: down-going P")
      pdf.savefig()
      plt.close()

      #-- plot SH radiation pattern
      Ash_max = max(np.max(np.abs(Ashu)), np.max(np.abs(Ashd)))
      
      # up-going SH
      fig = plt.figure()
      # contour
      ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect='equal')
      levels = np.linspace(-1.0,1.0,21)
      cs = ax.contour(xx, yy, Ashu/Ash_max, levels, colors='k')
      ax.clabel(cs, fontsize='x-small', inline=1, fmt='%.1f')
      # axis
      ax.plot(x_peri, y_peri, 'k')
      for theta in np.arange(10, 90, 10):
        xc = theta * np.sin(np.deg2rad(az))
        yc = theta * np.cos(np.deg2rad(az))
        ax.plot(xc, yc, '-', color=(0.5,0.5,0.5), linewidth=0.5)
      for theta in np.arange(0, 360, 10):
        xc = 90.0 * np.sin(np.deg2rad(theta))
        yc = 90.0 * np.cos(np.deg2rad(theta))
        ax.plot([0, xc], [0, yc], '-', color=(0.5,0.5,0.5), linewidth=0.5)
      # plot stations
      for station_id in station_dict:
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if dist_degree > min(dist) and dist_degree < max(dist):
          for arr in meta['ttime']:
            if arr.name[0] == "s":
              takeoff_angle = arr.takeoff_angle
              inc = 180.0 - takeoff_angle
              xc = inc * np.sin(np.deg2rad(azimuth))
              yc = inc * np.cos(np.deg2rad(azimuth))
              ax.plot(xc, yc, 'r+')
      # set axis
      ax.set_title("Radiation pattern: up-going SH")
      pdf.savefig()
      plt.close()

      #-- down-going SH
      Ash_max = max(np.max(np.abs(Ashu)), np.max(np.abs(Ashd)))
      fig = plt.figure()
      # contour
      ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect='equal')
      levels = np.linspace(-1.0,1.0,21)
      cs = ax.contour(xx, yy, Ashd/Ash_max, levels, colors='k')
      ax.clabel(cs, fontsize='x-small', inline=1, fmt='%.1f')
      # axis
      ax.plot(x_peri, y_peri, 'k')
      for theta in np.arange(10, 90, 10):
        xc = theta * np.sin(np.deg2rad(az))
        yc = theta * np.cos(np.deg2rad(az))
        ax.plot(xc, yc, '-', color=(0.5,0.5,0.5), linewidth=0.5)
      for theta in np.arange(0, 360, 10):
        xc = 90.0 * np.sin(np.deg2rad(theta))
        yc = 90.0 * np.cos(np.deg2rad(theta))
        ax.plot([0, xc], [0, yc], '-', color=(0.5,0.5,0.5), linewidth=0.5)
      # plot stations
      for station_id in station_dict:
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if dist_degree > min(dist) and dist_degree < max(dist):
          for arr in meta['ttime']:
            if arr.name[0] == "S":
              takeoff_angle = arr.takeoff_angle
              inc = takeoff_angle
              xc = inc * np.sin(np.deg2rad(azimuth))
              yc = inc * np.cos(np.deg2rad(azimuth))
              ax.plot(xc, yc, 'r+')
      # set axis
      ax.set_title("Radiation pattern: down-going SH")
      pdf.savefig()
      plt.close()

      #-- plot SV radiation pattern
      Asv_max = max(np.max(np.abs(Asvu)), np.max(np.abs(Asvd)))
      # up-going SV
      fig = plt.figure()
      ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect='equal')
      levels = np.linspace(-1.0,1.0,21)
      cs = ax.contour(xx, yy, Asvu/Asv_max, levels, colors='k')
      ax.clabel(cs, fontsize='x-small', inline=1, fmt='%.1f')
      ax.plot(x_peri, y_peri, 'k')
      for theta in np.arange(10, 90, 10):
        xc = theta * np.sin(np.deg2rad(az))
        yc = theta * np.cos(np.deg2rad(az))
        ax.plot(xc, yc, '-', color=(0.5,0.5,0.5), linewidth=0.5)
      for theta in np.arange(0, 360, 10):
        xc = 90.0 * np.sin(np.deg2rad(theta))
        yc = 90.0 * np.cos(np.deg2rad(theta))
        ax.plot([0, xc], [0, yc], '-', color=(0.5,0.5,0.5), linewidth=0.5)
      # plot stations
      for station_id in station_dict:
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if dist_degree > min(dist) and dist_degree < max(dist):
          for arr in meta['ttime']:
            if arr.name[0] == "s":
              takeoff_angle = arr.takeoff_angle
              inc = 180.0 - takeoff_angle
              xc = inc * np.sin(np.deg2rad(azimuth))
              yc = inc * np.cos(np.deg2rad(azimuth))
              ax.plot(xc, yc, 'r+')
      # set axis
      ax.set_title("Radiation pattern: up-going SV")
      pdf.savefig()
      plt.close()

      # down-going SV
      fig = plt.figure()
      ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect='equal')
      levels = np.linspace(-1.0,1.0,21)
      cs = ax.contour(xx, yy, Asvd/Asv_max, levels, colors='k')
      ax.clabel(cs, fontsize='x-small', inline=1, fmt='%.1f')
      ax.plot(x_peri, y_peri, 'k')
      for theta in np.arange(10, 90, 10):
        xc = theta * np.sin(np.deg2rad(az))
        yc = theta * np.cos(np.deg2rad(az))
        ax.plot(xc, yc, '-', color=(0.5,0.5,0.5), linewidth=0.5)
      for theta in np.arange(0, 360, 10):
        xc = 90.0 * np.sin(np.deg2rad(theta))
        yc = 90.0 * np.cos(np.deg2rad(theta))
        ax.plot([0, xc], [0, yc], '-', color=(0.5,0.5,0.5), linewidth=0.5)
      # plot stations
      for station_id in station_dict:
        station = station_dict[station_id]
        if station['stat']['code'] < 0:
          continue
        meta = station['meta']
        azimuth = meta['azimuth']
        dist_degree = meta['dist_degree']
        if dist_degree > min(dist) and dist_degree < max(dist):
          for arr in meta['ttime']:
            if arr.name[0] == "S":
              takeoff_angle = arr.takeoff_angle
              inc = takeoff_angle
              xc = inc * np.sin(np.deg2rad(azimuth))
              yc = inc * np.cos(np.deg2rad(azimuth))
              ax.plot(xc, yc, 'r+')
      # set axis
      ax.set_title("Radiation pattern: down-going SV")
      pdf.savefig()
      plt.close()


#
#======================================================
#

  def read_obs_grf(self,
      obs_dir='obs',
      syn_dir='syn', syn_band_code='MX', syn_suffix='.sem.sac',
      left_pad=100, right_pad=0, obs_preevent=100):
    """ read in observed seismograms and synthetic Green's functions.
    Note:
      1) left_pad: time length to pad before synthetics 
        right_pad: time length to pad after synthetics 
        obs_pretime: pre-event time length of obs (for noise assessment)
      2) use delta STF in simulation to approximate green's function
    """
    syn_orientation_codes = ['E', 'N', 'Z']

    event = self.data['event']
    station_dict = self.data['station']

    if left_pad < 0:
      warnings.warn("[arg] left_pad must g.e. 0, set to 0.0")
      left_pad = 0
    if right_pad < 0:
      right_pad = 0
      warnings.warn("[arg] right_pad must g.e. 0, set to 0.0")
    if obs_preevent < 0:
      obs_preevent = 50.0
      warnings.warn("[arg] obs_preevent must g.e. 0, set to 50")

    for station_id in station_dict:
      station = station_dict[station_id]
      meta = station['meta']
      channel = meta['channel']
 
      #------ get file paths of obs, syn seismograms
      obs_files = [ '{:s}/{:s}.{:s}'.format(
        obs_dir, station_id, x['code']) for x in channel ]
      syn_files = [ '{:s}/{:s}.{:2s}{:1s}{:s}'.format(
        syn_dir, station_id, syn_band_code, x, syn_suffix)
        for x in syn_orientation_codes ]

      #------ read in obs, syn seismograms
      try:
        obs_st  = read(obs_files[0])
        obs_st += read(obs_files[1])
        obs_st += read(obs_files[2])
        syn_st  = read(syn_files[0])
        syn_st += read(syn_files[1])
        syn_st += read(syn_files[2])
      except Exception as e:
        warn_str = "failed to read in sac files for %s, SKIP (%s)" \
            % (station_id, sys.exc_info()[0])
        warnings.warn(warn_str, UserWarning)
        print(e)
        station['stat']['code'] = -1
        station['stat']['msg'] = "failed to read in sac files [%s]" \
            % UTCDateTime.now().isoformat()
        continue

      #------ get time samples of syn seismograms
      if not is_equal( [ (tr.stats.starttime, tr.stats.delta, tr.stats.npts) \
          for tr in syn_st ] ):
        raise Exception('%s not equal time samples in'\
            ' synthetic seismograms.' % (station_id))
      tr = syn_st[0]
      syn_delta = tr.stats.delta
      syn_npts = tr.stats.npts
      # padding
      nl = int(left_pad/syn_delta)
      nr = int(right_pad/syn_delta)
      nt = syn_npts + nl + nr
      # ENZ_syn
      syn_starttime = tr.stats.starttime - nl * syn_delta
      syn_endtime = syn_starttime + (nt-1)*syn_delta
      syn_times = np.arange(nt) * syn_delta
      syn_ENZ = np.zeros((3, nt))
      for i in range(3):
        syn_ENZ[i,nl:(nl+syn_npts)] = syn_st[i].data

      #------ interpolate obs into the same time samples of syn
      obs_ENZ = np.zeros((3, nt))
      syn_nyq = 0.5/syn_delta
      flag = True
      for i in range(3):
        tr = obs_st[i]
        obs_npts = tr.stats.npts
        obs_delta = tr.stats.delta
        obs_starttime = tr.stats.starttime
        obs_endtime = tr.stats.endtime
        # check if obs has enough pre-event records
        first_arrtime = event['t0'] + meta['ttime'][0].time
        obs_tmin = first_arrtime - obs_preevent
        if obs_starttime > obs_tmin or obs_endtime < syn_endtime:
          flag = False
          warn_str = "%s: record not long enough, SKIP %s" \
              % (obs_files[i], station_id)
          warnings.warn(warn_str)
          print(obs_starttime, obs_endtime)
          print(obs_tmin, syn_endtime)
          station['stat']['code'] = -1
          station['stat']['msg'] = "%s [%s]" \
              % (warn_str, UTCDateTime.now().isoformat())
          break

        #DEBUG
        #print(station_id, obs_files[i]
        #obs_times = (obs_starttime-syn_starttime) + np.arange(obs_npts)*obs_delta
        #plt.plot(obs_times, tr.data, 'k')

        # lowpass below the nyquist frequency of synthetics
        # repeat twice to avoid numerical inaccuries
        tr.detrend(type='linear')
        tr.detrend(type='linear')
        # repeat lowpass filter twice to make sharper edge
        tr.filter('lowpass', freq=0.8*syn_nyq, corners=10, zerophase=True)
        tr.filter('lowpass', freq=0.8*syn_nyq, corners=10, zerophase=True)
        # interpolation: windowed sinc reconstruction
        obs_ENZ[i,:] = lanczos_interp1(tr.data, obs_delta,
            syn_times+(syn_starttime-obs_starttime), na=20)

        #DEBUG
        #plt.plot(obs_times, tr.data, 'r')
        ##plt.plot(syn_times, obs_ENZ[i,:], 'r')
        #plt.show()

      # if bad data, skip this station
      if not flag:
        continue

      #------ rotate obs to ENZ
      # projection matrix: obs = proj * ENZ => ENZ = inv(proj) * obs
      proj_matrix = np.zeros((3, 3))
      for i in range(3):
        chan = channel[i]
        sin_az = np.sin(np.deg2rad(chan['azimuth']))
        cos_az = np.cos(np.deg2rad(chan['azimuth']))
        sin_dip = np.sin(np.deg2rad(chan['dip']))
        cos_dip = np.cos(np.deg2rad(chan['dip']))
        # column vector = obs channel polarization 
        proj_matrix[i,0] = cos_dip*sin_az # proj to E
        proj_matrix[i,1] = cos_dip*cos_az # proj to N
        proj_matrix[i,2] = -sin_dip     # proj to Z
      # inverse projection matrix: ENZ = inv(proj) * obs
      inv_proj = np.linalg.inv(proj_matrix)
      obs_ENZ = np.dot(inv_proj, obs_ENZ)

      #------ record data 
      if 'waveform' not in station:
        station['waveform'] = {}
      waveform = station['waveform']
      waveform['time_sample'] = {
          'starttime': syn_starttime, 'delta': syn_delta,
          'nt': nt, 'nl': nl, 'nr': nr }
      waveform['obs'] = obs_ENZ
      waveform['grf'] = syn_ENZ

      station['stat']['code'] = 0
      station['stat']['msg'] = "read_obs_grf OK [%s]" \
          % (UTCDateTime.now().isoformat())

    #endfor station_id in station_dict:
  #enddef read_obs_grf

#
#======================================================
#

  def setup_window(self,
      window_list=[('F','p,P',[-30,50]), ('F','s,S',[-40,70])],
      filter_param=('butter', 3, [0.01, 0.10]),
      taper_param=('cosine', 0.1)):
    """ Setup data windows based on ray arrivals in 1D earth model.
    
      window_list: define data window
      [ (component, phases, [begin, end]), ...]

      filter_param: (type, freqlims)
    """
    # filter/taper parameters
    filter_dict = {'type': filter_param[0], 
        'order': filter_param[1], 'freqlim': filter_param[2]}

    if not 0.0 < taper_param[1] < 0.5:
      raise ValueError("taper ratio must lie between 0 and 0.5.")

    event = self.data['event']
    station_dict = self.data['station']

    # loop each station
    for station_id in station_dict:
      station = station_dict[station_id]
      meta = station['meta']
      arrivals = meta['ttime']
      baz = meta['back_azimuth']

      # initialize window dict
      station['window'] = {}
      window = station['window']

      # loop each window
      for win in window_list:
        comp = win[0]
        phase = win[1]
        signal_begin = float(win[2][0])
        signal_end = float(win[2][1])
        window_id = "%s.%s" % (comp, phase)

        # window time range
        phase_list = phase.split(',')
        ref_time = event['t0']
        ttime = []
        for arr in arrivals:
          if arr.name in phase_list:
            ttime.append(arr.time)
        if ttime:
          ref_time += min(ttime)
          #print("[INFO] phase %s: min(ttime)=%f, ref_time=%s" \
          #    % (phase, min(ttime), ref_time)
        else:
          print("[INFO] phase %s not found, use event origin time=%s" \
              % (phase, ref_time))
        starttime = ref_time + signal_begin
        endtime = ref_time + signal_end
        taper_dict = {'type':taper_param[0], 'ratio':taper_param[1],
            'starttime':starttime, 'endtime':endtime}

        # window polarity 
        if comp == 'Z': # vertcal component
          cmpaz = 0.0 
          cmpdip = -90.0
        elif comp == 'R': # radial component
          cmpaz = (baz + 180.0)%360.0
          cmpdip = 0.0
        elif comp == 'T': # tangential component (TRZ: right-hand convention)
          cmpaz = (baz - 90.0)%360.0
          cmpdip = 0.0
        elif comp == 'H': # horizontal particle motion 
          cmpaz = float('nan')
          cmpdip = 0.0
        elif comp == 'F': # 3-d particle motion 
          cmpaz = float('nan')
          cmpdip = float('nan')
        else:
          print("[WARN] %s: unrecognized component, SKIP." % (comp))
          continue
        polarity_dict = {'component':comp, 'azimuth': cmpaz, 'dip': cmpdip }

        # add window
        window[window_id] = {
          'stat': {
            'code': 0, 
            'msg': "created on "+UTCDateTime.now().isoformat() },
          'filter': filter_dict,
          'taper': taper_dict,
          'polarity': polarity_dict }

      #endfor win in window_list:
    #endfor station_id, station in station_dict.iteritems():

  #enddef setup_windows

#
#======================================================
#

  def measure_adj(self,
      plot=False,
      cc_delta=0.01, 
      weight_param={
        'SNR':[10,15], 
        'CC0':[0.5,0.7], 
        'CCmax':None,
        'dist':None}
      ):
    """ calculate adjoint sources (dchi_du)
        chi: misfit functional (normalized zero-lag correlation coef.)
        u: synthetic waveform
        weight_param: SNR, CCmax, CC0
    """
    #------
    event = self.data['event']
    station_dict = self.data['station']

    # loop each station
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue

      meta = station['meta']
      window_dict = station['window']

      waveform = station['waveform']
      time_sample = waveform['time_sample']
      syn_starttime = time_sample['starttime']
      syn_delta = time_sample['delta']
      syn_nyq = 0.5/syn_delta
      syn_nt = time_sample['nt']
      syn_nl = time_sample['nl']
      syn_nr = time_sample['nr']
      syn_times = syn_delta * np.arange(syn_nt)

      obs = waveform['obs']
      grf = waveform['grf']

      # source spectrum (moment-rate function)
      syn_freq = np.fft.rfftfreq(syn_nt, d=syn_delta)
      F_src = stf_gauss_spectrum(syn_freq, event['tau'])

      #------ loop each window
      # sum of adjoint sources from all windows
      dchi_du = np.zeros((3, syn_nt))
      dchi_dg = np.zeros((3, syn_nt))

      for window_id in window_dict:
        # window parameters
        window = window_dict[window_id]
        # skip bad windows
        if window['stat']['code'] < 0:
          continue

        #------ filter
        filter_dict = window['filter']
        filter_type = filter_dict['type']
        filter_order = filter_dict['order']
        filter_freqlim = filter_dict['freqlim']
        # filter design 
        if filter_type == 'butter':
          filter_b, filter_a = signal.butter(filter_order,
            np.array(filter_freqlim)/syn_nyq, btype='band')
        else:
          raise Exception("filter_type %s not recognized" % filter_type)
        # record filter coeff.
        filter_dict['a'] = filter_a
        filter_dict['b'] = filter_b

        #------ taper
        taper = window['taper']
        taper_type = taper['type']
        taper_ratio = taper['ratio']
        # time range
        window_starttime = taper['starttime']
        window_endtime = taper['endtime']
        window_len = window_endtime - window_starttime
        # taper design
        win_b = window_starttime - syn_starttime
        win_e = window_endtime - syn_starttime
        taper_width = window_len * min(taper_ratio, 0.5)
        win_c = [win_b, win_b+taper_width, win_e-taper_width, win_e]
        if taper_type == "cosine":
          win_func = cosine_taper(syn_times, win_c)
        else:
          raise Exception("taper_type not recognized.")
        taper['win'] = win_func

        #------ polarity 
        polarity = window['polarity']
        comp = polarity['component']
        cmpaz = polarity['azimuth']
        cmpdip = polarity['dip']
        if comp in ['Z', 'R', 'T']:
          sin_az = np.sin(np.deg2rad(cmpaz))
          cos_az = np.cos(np.deg2rad(cmpaz))
          sin_dip = np.sin(np.deg2rad(cmpdip))
          cos_dip = np.cos(np.deg2rad(cmpdip))
          n = np.array([ [cos_dip * sin_az], # cos(E, comp)
                 [cos_dip * cos_az], # N, comp
                 [-sin_dip] ])     # Z, comp
          proj_matrix = np.dot(n, n.transpose())
        elif comp == 'H': # horizontal vector 2d
          proj_matrix = np.identity(3)
          proj_matrix[2,2] = 0.0 # zero Z component
        elif comp == 'F': # full 3d vector
          proj_matrix = np.identity(3)
        else:
          print('[WARNING] %s:%s unrecognized component code, SKIP' \
              % (station_id, window_id))
          continue
        polarity['proj_matrix'] = proj_matrix

        #------ filter obs, syn
        #NOTE: use lfilter (causal filter) to avoid contamination from the right
        # end of the signal, but with asymmetric response and 
        # peak shift ~ 1/4 min. period (e.g. 0.01-0.1Hz -> 2.5s peak shift)
        # , however the duration of the filter response is determined by the
        # max. period (e.g. 0.01-0.1Hz -> ~50s). So the time window chosen 
        # should not be affected by the relatively small peak shift.
        #-- F * d
        obs_filt = signal.filtfilt(filter_b, filter_a, obs)
        #-- F * u (u = S*grf)
        syn_filt = signal.filtfilt(filter_b, filter_a, grf)
        syn_filt = np.fft.irfft(F_src*np.fft.rfft(syn_filt), syn_nt)
        #DEBUG
        #for i in range(3):
        #  plt.subplot(311+i)
        #  plt.plot(syn_times, obs[i,:], 'k', syn_times, obs_filt[i,:], 'r')
        #plt.show()
        #-- noise: use signals 40s before first arrival time on obs
        first_arrtime = event['t0'] + meta['ttime'][0].time
        #FIXME: better choice of the time length before first arrival? 
        tnoise = (first_arrtime - 40.0) - syn_starttime
        noise_idx = syn_times < tnoise
        #t = syn_times[noise_idx]
        #b = t[0]
        #e = t[-1]
        #taper_width = (e-b) * 0.1
        #win_c = [b, b+taper_width, e-taper_width, e]
        #taper = cosine_taper(t, win_c)
        # F * noise
        noise_filt = obs_filt[:,noise_idx]

        #------ apply window taper and polarity projection
        # obs = w * F * d
        obs_filt_win = np.dot(proj_matrix, obs_filt) * win_func 
        # syn = w * F * u (u = S*g)
        syn_filt_win = np.dot(proj_matrix, syn_filt) * win_func 
        # noise (only projection)
        noise_filt_win = np.dot(proj_matrix, noise_filt)
        #DEBUG
        #diff = obs_ENZ_win - syn_ENZ_win
        #for i in range(3):
        #  plt.subplot(311+i)
        #  plt.plot(syn_times, obs_ENZ_win[i,:], 'k')
        #  plt.plot(syn_times, syn_ENZ_win[i,:], 'r')
        #  plt.plot(syn_times, diff[i,:], 'c')
        #plt.show()

        #------ measure SNR (based on maximum amplitude)
        Amax_obs = np.sqrt(np.max(np.sum(obs_filt_win**2, axis=0)))
        Amax_syn = np.sqrt(np.max(np.sum(syn_filt_win**2, axis=0)))
        Amax_noise =  np.sqrt(np.max(np.sum(noise_filt_win**2, axis=0)))
        if Amax_obs == 0: # bad record
          print('[WARN] %s:%s:%s empty obs trace, SKIP.' \
              % (event_id, station_id, window_id))
          window['stat']['code'] = -1
          window['stat']['msg'] = "Amax_obs=0"
          continue
        if Amax_noise == 0: # could occure when the data begin time is too close to the first arrival
          print('[WARN] %s:%s:%s empty noise trace, SKIP.' \
              % (event_id, station_id, window_id))
          window['stat']['code'] = -1
          window['stat']['msg'] = "Amax_noise=0"
          continue
        snr = 20.0 * np.log10(Amax_obs/Amax_noise)
 
        #------ measure CC time shift (between w*F*d and w*F*u)
        obs_norm = np.sqrt(np.sum(obs_filt_win**2))
        syn_norm = np.sqrt(np.sum(syn_filt_win**2))
        # window normalization factor (without dt)
        Nw = obs_norm * syn_norm
        # NOTE the order (obs,syn) is important. The positive time on 
        # CC means shifting syn in the positive time direction to match
        # the observed obs, and vice verser.
        # [-(nt-1), nt) * dt
        cc = np.zeros(2*syn_nt-1)
        for i in range(3):
          cc += signal.fftconvolve(
              obs_filt_win[i,:], syn_filt_win[i,::-1], 'full')
        cc /= Nw
        #DEBUG
        #print(window_id
        #print(cc[syn_nt-2] - np.sum(obs_ENZ_win * syn_ENZ_win)
        #print(cc[syn_nt-1] - np.sum(obs_ENZ_win * syn_ENZ_win)
        #print(cc[syn_nt] - np.sum(obs_ENZ_win * syn_ENZ_win)
        #print(cc[syn_nt+1] - np.sum(obs_ENZ_win * syn_ENZ_win)
        #print(cc[syn_nt+2] - np.sum(obs_ENZ_win * syn_ENZ_win)
        #-- zero-lag cc coeff.
        CC0 = cc[syn_nt-1] #the n-th point corresponds to zero lag time 
        AR0 = CC0 * syn_norm / obs_norm # amplitude ratio syn/obs 
        #DEBUG
        #print(CC0 - np.sum(obs_ENZ_win * syn_ENZ_win)/obs_norm/syn_norm
        #-- interpolate cc to finer time samples
        CC_shift_range = window_len/2.0 #TODO: more reasonable choice?
        ncc = int(CC_shift_range / cc_delta)
        cc_times = np.arange(-ncc,ncc+1) * cc_delta
        if syn_delta < cc_delta:
          warnings.warn("syn_delta(%f) < cc_time_step(%f)" \
              % (syn_delta, cc_delta))
        ti = cc_times + (syn_nt-1)*syn_delta  # -(npts-1)*dt: begin time in cc
        cci = lanczos_interp1(cc, syn_delta, ti, na=20)
        # time shift at the maximum correlation
        imax = np.argmax(cci)
        CC_time_shift = cc_times[imax]
        CCmax = cci[imax]
        ARmax = CCmax * syn_norm / obs_norm # amplitude ratio: syn/obs

        #------ window weighting based on SNR and misfit
        weight = 1.0
        if 'SNR' in weight_param:
          weight *= cosine_taper(snr, weight_param['SNR'])
        if 'CCmax' in weight_param:
          weight *= cosine_taper(CCmax, weight_param['CCmax'])
        if 'CC0' in weight_param:
          weight *= cosine_taper(CC0, weight_param['CC0'])
        if 'dist' in weight_param:
          dist = meta['dist_degree']
          dist_range = np.array(weight_param['dist'])
          weight *= cosine_taper(dist, dist_range)

        #------ measure adjoint source
        # adjoint source: dchiw_du (misfit functional: zero-lag cc coef.)
        # dchiw_du = conj(F * [S]) * w * [ w * F * d - A * w * F * S * g] / N, 
        # , where A = CC0(un-normalized) / norm(u)**2, N = norm(d)*norm(u)
        Aw = CC0 * obs_norm / syn_norm # window amplitude raito
        #-- dchiw_du
        #NOTE: *dt is put back to Nw
        dchiw_du1 = win_func * (obs_filt_win - Aw*syn_filt_win) / Nw / syn_delta
        # apply conj(F), equivalent to conj(F*conj(adj))
        # for two-pass filter (zero phase) conj(F) = F
        dchiw_du = signal.filtfilt(filter_b, filter_a, dchiw_du1[:,::-1])
        dchiw_du = dchiw_du[:,::-1]
        #DEBUG
        #for i in range(3):
        #  plt.subplot(311+i)
        #  plt.plot(syn_times, dchiw_du1[i,:], 'k')
        #  plt.plot(syn_times, dchiw_du[i,:], 'r')
        #plt.show()
        # add into total dchi_du
        dchi_du += weight * dchiw_du
        #-- dchiw_dg = conj(S) * dchiw_du
        dchiw_dg = np.fft.irfft(np.conjugate(F_src) * np.fft.rfft(dchiw_du), syn_nt)
        # add into total dchi_dg
        dchi_dg += weight * dchiw_dg 
        #DEBUG
        #for i in range(3):
        #  plt.subplot(311+i)
        #  plt.plot(syn_times, dchiw_du[i,:], 'k')
        #  plt.plot(syn_times, dchiw_dg[i,:], 'r')
        #plt.show()

        #------ record results
        quality_dict = {
            'Amax_obs': Amax_obs, 'Amax_syn': Amax_syn, 
            'Amax_noise': Amax_noise, 'SNR': snr}
        cc_dict = {
            'time': cc_times, 'cc': cci,
            'cc_tshift': CC_time_shift,
            'CC0': CC0, 'CCmax': CCmax,
            'AR0': AR0, 'ARmax': ARmax,
            'Nw':Nw, 'Aw':Aw }
        window['quality'] = quality_dict
        window['cc'] = cc_dict
        window['weight'] = weight
        window['stat'] = {'code': 1, 
            'msg': "measure adj on "+UTCDateTime.now().isoformat()}

        #------ plot measure window and results 
        if plot:
          syn_npts = syn_nt - syn_nl - syn_nr
          syn_orientation_codes = ['E', 'N', 'Z']
          adj = dchiw_dg
          Amax_adj = np.sqrt(np.max(np.sum(adj**2, axis=0)))
          t = syn_times
          for i in range(3):
            plt.subplot(411+i)
            if i == 0:
              plt.title('%s dt %.2f CCmax %.3f ARmax %.3f CC0 %.3f '
                  'AR0 %.3f \nAobs %g Anoise %g SNR %.1f weight %.3f'
                  % (station_id, CC_time_shift, CCmax, ARmax, 
                    CC0, AR0, Amax_obs, Amax_noise, snr, weight) )
            idx_plt = range(syn_nl,(syn_nl+syn_npts))
            plt.plot(t[idx_plt], obs_filt[i,idx_plt]/Amax_obs, 'k', linewidth=0.2)
            plt.plot(t[idx_plt], syn_filt[i,idx_plt]/Amax_obs*Aw, 'r', linewidth=0.2)
            #plt.plot(t[noise_idx], noise_filt[i,:]/Amax_obs, 'b', linewidth=1.0)
            idx = (win_b <= syn_times) & (syn_times <= win_e)
            plt.plot(t[idx], obs_filt_win[i,idx]/Amax_obs, 'k', linewidth=1.0)
            plt.plot(t[idx], syn_filt_win[i,idx]/Amax_obs*Aw, 'r', linewidth=1.0)
            plt.plot(t[idx_plt], adj[i,idx_plt]/Amax_adj, 'c', linewidth=1.0)
            plt.ylim((-1.5, 1.5))
            #plt.xlim((min(t), max(t)))
            plt.xlim((t[syn_nl], t[syn_nl+syn_npts-1]))
            plt.ylabel(syn_orientation_codes[i])
          plt.subplot(414)
          plt.plot(cc_times, cci, 'k-')
          plt.xlim((min(cc_times), max(cc_times)))
          plt.ylabel(window_id)
          plt.show()
      #====== end for window_id in windows:

      #------ store adjoint source for this station
      station['dchi_du'] = dchi_du
      station['dchi_dg'] = dchi_dg

      #DEBUG
      #for i in range(3):
      #  plt.subplot(311+i)
      #  plt.plot(syn_times, dchi_du[i,:], 'k')
      #  plt.plot(syn_times, dchi_dg[i,:], 'r')
      #plt.show()

    #endfor station_id in station_dict:
  #enddef measure_windows_for_one_station(self,

#
#======================================================
#

  def output_adj(self, 
      adj_type='dchi_dg',
      out_dir='adj',
      syn_band_code='MX'):
    """Output adjoint sources
    NOTE:
      1) dchi_dg: use tau=0 in forward/adjoint simulation
      1) dchi_du: use real tau in forward/adjoint simulation
    """
    syn_orientation_codes = ['E', 'N', 'Z']
    event = self.data['event']
    station_dict = self.data['station']

    tr = Trace()
    for station_id in station_dict:
      station = station_dict[station_id]
      # skip rejected statations
      if station['stat']['code'] < 0:
        continue
      # time samples
      waveform = station['waveform']
      time_sample = waveform['time_sample']
      syn_starttime = time_sample['starttime']
      syn_delta = time_sample['delta']
      syn_nt = time_sample['nt']
      syn_nl = time_sample['nl']
      syn_nr = time_sample['nr']
      # without padding
      npts = syn_nt - syn_nl - syn_nr
      starttime = syn_starttime + syn_nl*syn_delta
      # time samples for ascii output, referred to origin time
      syn_times = np.arange(npts)*syn_delta
      b = starttime - event['t0']
      syn_times += b

      # adjoint source
      if adj_type == 'dchi_du':
        adj = station['dchi_du']
      elif adj_type == 'dchi_dg':
        adj = station['dchi_dg']
      else:
        raise Exception('unknown adj_type: %s (dchi_du or dchi_dg) ' \
            % (adj_type))

      # loop ENZ
      for i in range(3):
        tr.data = adj[i, syn_nl:(syn_nl+npts)]
        tr.stats.starttime = starttime
        tr.stats.delta = syn_delta

        out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
            out_dir, station_id, syn_band_code,
            syn_orientation_codes[i])

        # sac format
        tr.write(out_file + '.adj.sac', 'sac')

        # ascii format (needed by SEM)
        # time is relative to event origin time: t0
        with open(out_file+'.adj','w') as fp:
          for j in range(npts):
            fp.write("{:16.9e}  {:16.9e}\n".format(
              syn_times[j], adj[i,syn_nl+j]))

      #endfor i in range(3):
    #endfor station_id in station_dict:
  #enddef output_adjoint_source

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

#  def waveform_der_stf(self):
#    """ Calculate waveform derivatives for source time function (t0, tau)
#    """
#    event = self.data['event']
#    t0 = event['t0']
#    tau = event['tau']
#
#    station_dict = self.data['station']
#    for station_id in station_dict:
#      station = station_dict[station_id]
#      # skip rejected statations
#      if station['stat']['code'] < 0:
#        continue
#
#      # source time function
#      waveform = station['waveform']
#      time_sample = waveform['time_sample']
#      starttime = time_sample['starttime']
#      dt = time_sample['delta']
#      nt = time_sample['nt']
#      t = np.arange(nt) * dt + (starttime - t0) #referred to t0
#      # s(t), Ds(t)/Dt0, Ds(t)/Dtau
#      freq = np.fft.rfftfreq(nt, d=dt)
#      F_src, F_ds_dt0, F_ds_dtau = stf_gauss_spectrum(freq, tau)
#
#      #------ waveform derivative
#      # green's function
#      grf = waveform['grf']
#      # convolve Ds(t)/Dt0,tau with Green's function
#      du_dt0 = np.fft.irfft(F_ds_dt0 * np.fft.rfft(grf), nt)
#      du_dtau = np.fft.irfft(F_ds_dtau * np.fft.rfft(grf), nt)
#      # zero records before origin time (wrap around from the end)
#      idx = t < -5.0*tau
#      du_dt0[:,idx] = 0.0
#      du_dtau[:,idx] = 0.0
#
#      #------ misfit derivative
#      # adjoint source = Dchi/Du
#      dchi_du = station['dchi_du']
#      dchi_dt0 = np.sum(dchi_du * du_dt0) * dt
#      dchi_dtau = np.sum(dchi_du * du_dtau) * dt
#
#      #------ record derivatives
#      if 'waveform_der' not in station:
#        station['waveform_der'] = {}
#      station['waveform_der']['dt0'] = {
#          'dm':1.0, 'du':du_dt0, 'dchi':dchi_dt0 }
#      station['waveform_der']['dtau'] = {
#          'dm':1.0, 'du': du_dtau, 'dchi':dchi_dtau }
#
#      # DEBUG
#      #print(dchi_dt0, dchi_dtau
#      #for i in range(3):
#      #  plt.subplot(311+i)
#      #  #plt.plot(t, dchi_du[i,:], 'k')
#      #  plt.plot(t, du_dt0[i,:], 'b', t, du_dtau[i,:], 'r')
#      #plt.show()
#
#  #enddef derivative_stf(self)

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
      raise Exception('src_frechet not set.')
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
      #print(dchi
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
      #print(dchi
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
        obs_filt = signal.filtfilt(filter_b, filter_a, obs)
        # F * u (u = S*grf)
        syn_filt = signal.filtfilt(filter_b, filter_a, grf)
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
        #print("Nw: %e %e" % (Nw, window['cc']['Nw'])
        #print("Aw: %e %e" % (Aw, window['cc']['Aw'])

        #------ filter differential seismograms (w * F * du_dm)
        wFdu = {}
        for param in src_param:
          du = waveform_der[param]['du']
          Fdu = signal.filtfilt(filter_b, filter_a, du)
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

#   print("dchi_dm:"
#   print(dchi_dm 

#   print("hessian:"
#   print(hessian

#   print("====== 0:4:"
#   w, v = np.linalg.eigh(hessian, UPLO='U')
#   print(w
#   print(v
#   x, residual, rank, sigval = np.linalg.lstsq(hessian, -dchi_dm)
#   print(" inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print("dt0: \n", x[0]
#   print("dtau:\n", x[1]
#   print("dxs: \n", x[2]*self.data['src_perturb']['xs'] 
#   print("dmt: \n", x[3]*self.data['src_perturb']['mt'] 

#   print("====== only 0:3"
#   h3 = hessian[0:3,0:3]
#   v3 = dchi_dm[0:3]
#   w, v = np.linalg.eigh(h3, UPLO='U')
#   print(w
#   print(v
#   x, residual, rank, sigval = np.linalg.lstsq(h3, -v3)
#   print("inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print("dt0: \n", x[0]
#   print("dtau:\n", x[1]
#   print("dxs: \n", x[2]*self.data['src_perturb']['xs'] 
#   #print("dmt: \n", x[3]*self.data['src_perturb']['dmt'] 

#   print("====== only 0:2"
#   h3 = hessian[0:2,0:2]
#   v3 = dchi_dm[0:2]
#   w, v = np.linalg.eigh(h3, UPLO='U')
#   print(w
#   print(v
#   x, residual, rank, sigval = np.linalg.lstsq(h3, -v3)
#   print("inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print("dt0: \n", x[0]
#   print("dtau:\n", x[1]

#   print("====== only 0,2"
#   idx = [0,2]
#   hess = hessian[idx,:][:,idx]
#   kernel = dchi_dm[idx]
#   print(hess
#   eigs, eigv = np.linalg.eigh(hess, UPLO='U')
#   print(eigs
#   print(eigv
#   x, residual, rank, sigval = np.linalg.lstsq(hess, -kernel)
#   print("inv(hessian)*(-1.0 * dchi_dm): \n", x
#   print("dt0: \n", x[0]
#   print("dxs:\n", x[1]

# #enddef measure_windows_for_one_station(self,

#
#======================================================
#

  def cc_perturbed_seismogram(self,
      dm={'t0':None, 'xs':None},
      dist_min=0.0,
      dist_max=180.0,
      plot=False
      ):
    """ calculate normalized zero-lag cc for perturbed seismograms from linear combination of waveform derivatives

    dm={<model_name>:<model_vector>, ...}
    model_name: ['dt0', 'dxs']
    model_vector: ndarray of the same length

    dist_range : array(2)
      if set, use only stations within the specified distance range (degree)

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
      # skip rejected stations
      if station['stat']['code'] < 0:
        continue

      # skip stations out of the selected distance range
      gcarc = station['meta']['dist_degree']
      if gcarc < dist_min or gcarc > dist_max:
        #print("SKIP %s: gcarc=%f not in (%f, %f)" % (
        #  station_id, gcarc, dist_min, dist_max))
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
        # skip window with zero weight
        if np.isclose(weight, 0.0):
          continue
        weight_sum += weight
        # filter
        filter_dict = window['filter']
        filter_a = filter_dict['a']
        filter_b = filter_dict['b']
        # taper
        win_func = window['taper']['win']
        win_starttime = window['taper']['starttime']
        win_endtime = window['taper']['endtime']
        if plot:
          syn_times = np.arange(syn_nt) * syn_delta
          plt.plot(syn_times, win_func)
          plt.show()
          plt.title("wind_func")
        # polarity projection 
        proj_matrix = window['polarity']['proj_matrix']
        #-- filter,project,taper obs
        # F * d
        obs_filt = signal.filtfilt(filter_b, filter_a, obs)
        # w * p * F * d (window,project,filter)
        wpFd = np.dot(proj_matrix, obs_filt) * win_func 
        norm_wpFd = np.sqrt(np.sum(wpFd**2))
        #-- filter,project grf 
        # F * g
        grf_filt = signal.filtfilt(filter_b, filter_a, grf)
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
            dg_filt = signal.filtfilt(filter_b, filter_a, dg)
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
            obs_filt_proj = np.dot(proj_matrix, obs_filt)
            syn_filt_proj = np.dot(proj_matrix, syn_filt)
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
              # whole trace, projected
              plt.plot(syn_times[idx_plt], obs_filt_proj[i,idx_plt]/Amax_obs, 
                  'k', linewidth=0.2)
              plt.plot(syn_times[idx_plt], syn_filt_proj[i,idx_plt]/Amax_syn,
                  'b', linewidth=0.2)
              # windowed trace
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

  def search1d_cc_perturbed_seismogram(self,
      dm_range = {
        't0': [-5,5], 
        'tau': [-5,5],
        'xs': [-5,5],
        'mt': [-5,5],
        },
      dist_min=0.0,
      dist_max=180.0,
      ngrid=4,
      range_ratio=0.7,
      plot_seism=False,
      ):
    """ 
    calculate misfit over 2D model grids based on perturbed seismograms 

    Parameters
    ----------
    range_ratio : scalar
      between (0, 1), used to shrink the search range to cc values above range_ratio
      between cc_max and cc_min.
    """
    # grid parameters
    model_name = [x for x in dm_range]

    range_ratio = float(range_ratio)
    if range_ratio < 0.5 or range_ratio > 0.9:
      raise ValueError("range_ratio should be within (0.5, 0.9)")

    dm_opt = {}
    for par in model_name:
      dm_opt[par] = 0.0

    niter = 0
    while niter < 10:
      print("====== iteration: %d" % (niter))

      # loop each parameter
      for par in model_name:
        print("------ %s" % (par))
        xlim = np.array(dm_range[par])

        if xlim.size == 1:
          dm_opt[par] = xlim[0]
          print("fixed at %.5e" % (dm_opt[par]))
          continue
        else:
          print("range: %.5e, %.5e" % (np.min(xlim), np.max(xlim)))
          xx = np.linspace(np.min(xlim), np.max(xlim), ngrid)

          # search
          dm_grid = {}
          dm_grid[par] = xx
          for par_fix in model_name:
            if par_fix != par:
              dm_grid[par_fix] = np.ones(ngrid) * dm_opt[par_fix]
          cc, weight = self.cc_perturbed_seismogram(
              dm=dm_grid, 
              dist_min=dist_min, 
              dist_max=dist_max,
              plot=plot_seism)
          cc /= weight # weighted average of normalized zero-lag CC

          # interpolation
          x_i = np.linspace(np.min(xlim), np.max(xlim), 100)
          interp = interpolate.interp1d(xx, cc, kind='cubic')
          cc_i = interp(x_i)
          imax = np.argmax(cc_i)
          x_max = x_i[imax]
          cc_max = cc_i[imax]

          #change search range
          imin = np.argmin(cc_i)
          cc_min = cc_i[imin]
          cc_thred = cc_min + (cc_max-cc_min)*range_ratio
          x1 = x_i[cc_i>cc_thred]
          dm_range[par] = [np.min(x1), np.max(x1)]

          dm_opt[par] = x_max
          print("maximum at %.5e with cc %.5e" % (dm_opt[par], cc_max))

      # print(out optimal model
      print("------ optimal alpha")
      for par in model_name:
        print("%s: %.5e" %(par, dm_opt[par]))
      niter += 1

      sys.stdout.flush()

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

    zz_min = np.min(zz)

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

          xx_1d = x[ix]
          x_i = np.linspace(np.min(xx_1d), np.max(xx_1d), 100)

          interp = interpolate.interp1d(xx_1d, zz_1d, kind='cubic')
          zz_i = interp(x_i)
          imax_i = np.argmax(zz_i)
          x_i_max = x_i[imax_i]
          zz_i_max = zz_i[imax_i]

          # plot cross-sections through maximum point
          #fig = plt.figure(figsize=(8.5,11))
          fig = plt.figure()
          plt.plot(xx_1d, zz_1d, 'ro', x_i, zz_i)

          title_str = "Weighted average CC  (%.2e, %.2e)" % (x_i_max, zz_i_max)
          plt.title(title_str)

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
          print(xx_2d)
          print(yy_2d)
          print(zz_2d)
          if ix > iy:
            zz_2d = np.transpose(zz_2d)
            xx_2d = np.transpose(xx_2d)
            yy_2d = np.transpose(yy_2d)

          # make interpolation grid
          x_i = np.linspace(np.min(xx_2d), np.max(xx_2d), 100)
          y_i = np.linspace(np.min(yy_2d), np.max(yy_2d), 100)
          xx_i, yy_i = np.meshgrid(x_i, y_i)

          interp = interpolate.interp2d(xx_2d, yy_2d, zz_2d, kind='cubic')
          zz_i = interp(x_i, y_i)
          imax_i = np.argmax(zz_i)
          idx_max_i = np.unravel_index(imax_i, xx_i.shape)
          zz_i_max = zz_i[idx_max_i]
          x_i_max = xx_i[idx_max_i]
          y_i_max = yy_i[idx_max_i]

          # plot 2D surface through the maximum
          #fig = plt.figure(figsize=(8.5, 11))
          fig = plt.figure()

          levels = zz_min + (zz_max-zz_min)*np.linspace(0.5, 1.0, 20)
          #cs = plt.contour(xx_2d, yy_2d, zz_2d, levels, colors='k')
          cs = plt.contour(xx_i, yy_i, zz_i, levels, colors='k')
          plt.clabel(cs, fontsize=9, inline=1)

          #x_max = xx[ix][idx_max]
          #y_max = xx[iy][idx_max]
          #text_str = "(%.1f,%.1f)" % (x_max, y_max)
          #plt.text(x_max, y_max, text_str,
          #    horizontalalignment="center", verticalalignment="top")
          #plt.plot(x_max, y_max, 'ko')

          #text_str = "%.2e,%.2e" % (x_i_max, y_i_max)
          #plt.text(x_i_max, y_i_max, text_str,
          #    horizontalalignment="center", verticalalignment="top")
          plt.plot(x_i_max, y_i_max, 'k+')

          title_str = "Weighted average CC (%.2e,%.2e)" % (x_i_max, y_i_max)
          plt.title(title_str)

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
      print("[ERROR] %s does NOT exist. Exit" % (event_id))
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
      print("[ERROR] %s does NOT exist. Exit" \
          % (event_id))
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
      twopass_filt=False,
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
    #------ selection parameters
    plot_time = plot_param['time']
    plot_azbin = plot_param['azbin']
    plot_rayp = plot_param['rayp']
    plot_window_id = plot_param['window_id']

    plot_SNR = np.array(plot_param['SNR'])
    plot_CC0 = np.array(plot_param['CC0'])
    plot_CCmax = np.array(plot_param['CCmax'])
    plot_dist = np.array(plot_param['dist'])

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

      print(azmin, azmax)

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
          if dist_degree < np.min(plot_dist) or dist_degree > np.max(plot_dist):
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
        if plot_SNR and quality['SNR']<np.min(plot_SNR):
          continue
        if plot_CC0 and cc['CC0']<np.min(plot_CC0):
          continue
        if plot_CCmax and cc['CCmax']<np.min(plot_CCmax):
          continue

        # get seismograms: syn/obs
        waveform = station['waveform']
        time_sample = waveform['time_sample']
        syn_starttime = time_sample['starttime']
        syn_npts = time_sample['nt']
        syn_delta = time_sample['delta']
        syn_nyq = 0.5/syn_delta
        obs = waveform['obs']
        grf = waveform['grf']
        # filter parameter
        filter_param = window['filter']
        filter_a = filter_param['a']
        filter_b = filter_param['b']
        # filter seismograms 
        obs = signal.filtfilt(filter_b, filter_a, obs)
        grf = signal.filtfilt(filter_b, filter_a, grf)
        # convolve stf on grf
        syn_freq = np.fft.rfftfreq(syn_npts, d=syn_delta)
        F_src = stf_gauss_spectrum(syn_freq, event['tau'])
        syn = np.fft.irfft(F_src*np.fft.rfft(grf), syn_npts)
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
        t0 = syn_starttime - event['t0']
        # time of samples referred to centroid time
        syn_times = syn_delta*np.arange(syn_npts) + t0
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