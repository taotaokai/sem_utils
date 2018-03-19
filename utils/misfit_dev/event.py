#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import warnings
#
import numpy as np
#
from obspy import UTCDateTime
#from obspy.taup import TauPyModel
#from obspy.imaging.beachball import Beach
#
import pyproj
#
#from matplotlib import colors, ticker, cm
#from matplotlib.backends.backend_pdf import PdfPages 
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap


class Event:
  """
  Event info

  Unit: kg,m,s

  Coordinate: ECEF cartesian

  Attributes
  ----------
  code # event ID
  header, # have centroid time in it
  latitude, longitude, depth, 
  centroid_time,
  gauss_width:, # gaussian width exp(-(t-t0)^2/tau^2)/tau/pi^0.5
  xs: [x, y, z], # source location in ECEF coordinate
  mt, # moment tensor in ECEF coord (N*m)
  mt_rtp, # in spherical coord

  Methods
  -------
  read_cmtsolution 
  write_cmtsolution 

  """

#
#======================================================
#

  def __init__(self):
    self.code = "None"
    self.header = [] 
    self.longitude = 0.0
    self.latitude = 0.0 
    self.depth_km = 0.0 
    self.centroid_time = UTCDateTime(2000,1,1)
    self.gauss_width = 0.0
    self.xs = np.zeros(3)
    self.mt = np.zeros((3,3))
    self.mt_rtp = np.zeros((3,3))

#
#======================================================
#

  def read_cmtsolution(self, cmt_file, isECEF=True):
    """
    Read in event info from CMTSOLUTION file
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

    if isECEF:
      tau = float(lines[3][1])
      x = float(lines[4][1])
      y = float(lines[5][1])
      z = float(lines[6][1])
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
    header[1] = str(t0.year)
    header[2] = str(t0.month)
    header[3] = str(t0.day)
    header[4] = str(t0.hour)
    header[5] = str(t0.minute)
    header[6] = str(t0.second + 1.0e-6*t0.microsecond)

    # moment tensor
    # isECEF=false: 1,2,3 -> r,theta,phi
    # isECEF=true:  1,2,3 -> x,y,z
    m11 = float( lines[7][1])
    m22 = float( lines[8][1])
    m33 = float( lines[9][1])
    m12 = float(lines[10][1])
    m13 = float(lines[11][1])
    m23 = float(lines[12][1])
    mt = np.array([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]])
    # transform from spherical to cartesian coordinate
    r = (x**2 + y**2 + z**2)**0.5
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    # rotation matrix
    sthe = np.sin(theta)
    cthe = np.cos(theta)
    sphi = np.sin(phi)
    cphi = np.cos(phi)
    # D(x,y,z)/D(r,t,p): rtp -> xyz 
    a = np.array(
        [ [ sthe*cphi, cthe*cphi, -1.0*sphi ],
          [ sthe*sphi, cthe*sphi,      cphi ],
          [ cthe     , -1.0*sthe,      0.0  ] ])
    if isECEF:
      mt_xyz = mt
      mt_rtp = np.dot(np.dot(np.transpose(a), mt), a)
    else: # spherical coordinate
      a = np.array(
          [ [ sthe*cphi, cthe*cphi, -1.0*sphi ],
            [ sthe*sphi, cthe*sphi,      cphi ],
            [ cthe     , -1.0*sthe,      0.0  ] ])
      # harvard cmt use dyn*cm, change to N*m
      mt_rtp = mt*1.0e-7
      mt_xyz = np.dot(np.dot(a, mt), np.transpose(a))

    # add event
    self.code = event_id
    self.header = header
    self.longitude = lon
    self.latitude = lat
    self.depth_km = dep
    self.centroid_time = t0
    self.gauss_width = tau
    self.xs = np.array([x, y, z])
    self.mt = mt_xyz
    self.mt_rtp = mt_rtp

#
#======================================================
#

  def write_cmtsolution(self, out_file, isECEF=True):
    """ 
    Output CMTSOLUTION file 
    """
    # modify origin time in header line to have centroid time 
    header = self.header
    t0 = self.centroid_time
    header[1] = str(t0.year)
    header[2] = str(t0.month)
    header[3] = str(t0.day)
    header[4] = str(t0.hour)
    header[5] = str(t0.minute)
    header[6] = str(t0.second + 1.0e-6*t0.microsecond)

    if isECEF:
      with open(out_file, 'w') as fp:
        fp.write('%s\n' % ' '.join(self.header))
        fp.write('%-18s %s\n' % ('event name:', self.code))
        fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
        fp.write('%-18s %+15.8E\n' % ('tau(s):',   self.gauss_width))
        fp.write('%-18s %+15.8E\n' % ('x(m):',     self.xs[0]))
        fp.write('%-18s %+15.8E\n' % ('y(m):',     self.xs[1]))
        fp.write('%-18s %+15.8E\n' % ('z(m):',     self.xs[2]))
        fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', self.mt[0,0]))
        fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', self.mt[1,1]))
        fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', self.mt[2,2]))
        fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', self.mt[0,1]))
        fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', self.mt[0,2]))
        fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', self.mt[1,2]))
    else:
      warnings.warn("isECEF=False not implemented now")