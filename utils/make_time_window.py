#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Generate time windows based purely on event/station locations
and correspoding phases. This serves as the first step in making
time windows. Later operations may apply to adjust the windows.


"""
import sys
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.core.util.geodetics import gps2DistAzimuth, kilometer2degrees

taup_model = TauPyModel(model="iasp91")

event_info = 'CMTSOLUTION'
station_info = 'IRIS.station'
winspec_list = 'winspec.list'

#====== read event info
with open(event_info, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]
    event = [x.split() for x in lines]

evyr  = event[0][1]
evmo  = event[0][2]
evday = event[0][3]
evhr  = event[0][4]
evmin = event[0][5]
evsec = event[0][6]
evdt = float(event[2][2])

isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(\
        evyr,evmo,evday,evhr,evmin,evsec)
evotime = UTCDateTime(isotime) + evdt 

evnm = event[1][2]
evla = float(event[4][1])
evlo = float(event[5][1])
evdp = float(event[6][1])

event = {'id':evnm, 'centroid time':evotime, 'latitude':evla,
'longitude':evlo, 'depth':evdp}

print '#event info: ', event

#====== read station metadata 
with open(station_info, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.replace('\n','').split('|') for x in lines]
stations = [ {'network':     x[0],
              'station':     x[1],
              'location':    x[2],
              'channel':     x[3],
              'latitude':    float(x[4]),
              'longitude':   float(x[5]),
              'elevation':   float(x[6]),
              'depth':       float(x[7]),
              'azimuth':     float(x[8]),
              'dip':         float(x[9]),
              'starttime':   UTCDateTime(x[15]),
              'endtime':     UTCDateTime(x[16]),
             } for x in lines ]

#====== read window specification list
with open(winspec_list, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.split() for x in lines]
winspecs = [ {'phase':   x[0].split(','),
              'channel': x[1],
              'begin':   float(x[2]),
              'end':     float(x[3]),
             } for x in lines ]

print '#window specification: ', winspecs

#====== make window selection list
for channel in stations:

    if (channel['channel'] != 'BHZ'): 
        continue

    if not(channel['starttime'] < evotime and channel['endtime'] > evotime):
        continue

    dist, az, baz = gps2DistAzimuth(event['latitude'], event['longitude'],
            channel['latitude'], channel['longitude'])

    for winspec in winspecs:

        arrivals = taup_model.get_travel_times(
                source_depth_in_km=event['depth'], 
                distance_in_degree=kilometer2degrees(dist/1000.0), 
                phase_list=winspec['phase'])

        if winspec['channel'] == 'BHZ':
            azimuth = 0.0 
            dip = -90.0
        elif winspec['channel'] == 'BHR':
            azimuth = (baz + 180.0)%360.0
            dip = 0.0
        elif winspec['channel'] == 'BHT':
            azimuth = (baz + 90.0)%360.0
            dip = 0.0
        else:
            azimuth = float('nan')
            dip = float('nan')

        if arrivals:
            print '{:s}|{:s}|{:s}|{:s}|{:.1f}|{:.1f}|{:s}|{:s}'.format( 
                    channel['network'], channel['station'],
                    channel['location'], winspec['channel'], azimuth, dip,
                    event['centroid time'] + arrivals[0].time + winspec['begin'],
                    event['centroid time'] + arrivals[0].time + winspec['end'])
        else:
            print "[ERROR] phase not found: ", winspec['phase']
