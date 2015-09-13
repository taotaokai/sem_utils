#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Create initial time window list

Generate time windows based purely on event/station locations
and correspoding phases. This serves as the first step in making
time windows. Later operations may apply to adjust the windows.

Window specification file:
<<<
#phase | orientation | begin | end
p,P F -30 50
p,P Z -30 50
p,P R -30 50
s,S F -50 70
s,S H -50 70
s,S Z -50 70
s,S R -50 70
s,S T -50 70
<<<

Output:
<<<
IC|MDJ|00|BHZ,BHE,BHN|0.0,90.0,0.0|-90.0,0.0,0.0|44.6170|129.5908|220.0|50.0|p,P|F|nan|nan|2010-02-18T01:14:06.665653Z|2010-02-18T01:15:26.665653Z
IC|MDJ|00|BHZ,BHE,BHN|0.0,90.0,0.0|-90.0,0.0,0.0|44.6170|129.5908|220.0|50.0|p,P|Z|0.0|-90.0|2010-02-18T01:14:06.665653Z|2010-02-18T01:15:26.665653Z
IC|MDJ|00|BHZ,BHE,BHN|0.0,90.0,0.0|-90.0,0.0,0.0|44.6170|129.5908|220.0|50.0|p,P|R|339.6|0.0|2010-02-18T01:14:06.665653Z|2010-02-18T01:15:26.665653Z
<<<

"""

import sys
import numpy as np
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.core.util.geodetics import gps2DistAzimuth, kilometer2degrees

#====== parameters
event_info = 'CMTSOLUTION'
station_file = 'IRIS.station.correctYP'
#station_info = 'CEA_BO_IRIS.metadata'
#stationxml_file = 'test.stationxml'
winspec_list = 'winspec.list'

output_file = 'IRIS.window'

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

event = {'code':evnm, 'centroid_time':evotime, 'latitude':evla,
'longitude':evlo, 'depth':evdp}

print '#event info: ', event

#====== read metadata at channel level
# read in stationxml
#station_inventory = read_inventory(stationxml_file, format='stationxml')
with open(station_file, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.replace('\n','').split('|') for x in lines]

stations = {}
for x in lines:
    net_sta_loc = (x[0], x[1], x[2])
    channel = {'channel':     x[3],
               'latitude':    float(x[4]),
               'longitude':   float(x[5]),
               'elevation':   float(x[6]),
               'depth':       float(x[7]),
               'azimuth':     float(x[8]),
               'dip':         float(x[9]),
               'starttime':   UTCDateTime(x[15]),
               'endtime':     UTCDateTime(x[16])}
    if net_sta_loc not in stations:
        stations[net_sta_loc] = []
    stations[net_sta_loc].append(channel)

#====== read window specification list
with open(winspec_list, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]

lines = [x.split() for x in lines]
winspecs = [ {'phase':       x[0],
              'orientation': x[1],
              'begin':   float(x[2]),
              'end':     float(x[3]),
             } for x in lines ]

print '#window specification: ', winspecs

#====== make window selection list
taup_model = TauPyModel(model="iasp91")

f = open(output_file, 'w')
print >> f, '#event info: ', event
print >> f, '#window specification: ', winspecs

for net_sta_loc in stations:

    # select channels which record the event 
    channels = [x for x in stations[net_sta_loc] if 
            x['starttime'] < evotime and x['endtime'] > evotime]

    # check number of components = 3
    if len(channels) != 3:
        print '[WARN] no exactly 3 components found for ', net_sta_loc
        continue

    # check same geo-coordinates fo the 3-components
    if not((channels[0]['latitude'] == \
            channels[1]['latitude'] == \
            channels[2]['latitude']) and \
           (channels[0]['longitude'] == \
            channels[1]['longitude'] == \
            channels[2]['longitude']) and \
           (channels[0]['elevation'] == \
            channels[1]['elevation'] == \
            channels[2]['elevation']) and \
           (channels[0]['depth'] == \
            channels[1]['depth'] == \
            channels[2]['depth'])):
        print '[WARN] not the same geo-coordinates for ', net_sta_loc
        continue

    # check channel orientations
    Z_comp = [(x['channel'], x['azimuth'], x['dip']) \
            for x in channels if x['channel'][2] == 'Z']

    if len(Z_comp) != 1 or abs(Z_comp[0][2]) != 90.0: 
        print '[WARN] Problem with Z channel: ', Z_comp, net_sta_loc
        continue

    H_comp = [(x['channel'], x['azimuth'], x['dip']) \
            for x in channels if x['channel'][2] != 'Z']

    if len(H_comp) != 2 or \
            abs(H_comp[0][2]) != 0.0 or abs(H_comp[1][2]) != 0.0 or \
            abs(np.cos(np.deg2rad(H_comp[0][1] - H_comp[1][1]))) > 0.1: 
        print '[WARN] Problem with horizontal channels: ', \
                H_comp, net_sta_loc, stations[net_sta_loc]
        continue

    #
    dist, az, baz = gps2DistAzimuth(event['latitude'], event['longitude'],
            channels[0]['latitude'], channels[0]['longitude'])

    for winspec in winspecs:

        orientation_code = winspec['orientation']
        if orientation_code == 'Z':
            winspec['azimuth'] = 0.0 
            winspec['dip'] = -90.0
        elif orientation_code == 'R':
            winspec['azimuth'] = (baz + 180.0)%360.0
            winspec['dip'] = 0.0
        elif orientation_code == 'T':
            winspec['azimuth'] = (baz + 90.0)%360.0
            winspec['dip'] = 0.0
        elif orientation_code == 'H':
            winspec['azimuth'] = float('nan')
            winspec['dip'] = 0.0
        elif orientation_code == 'F':
            winspec['azimuth'] = float('nan')
            winspec['dip'] = float('nan')
        else:
            print '[ERROR] unrecognized orientation code: ', \
                    winspec['orientation']
            sys.exit()

        arrivals = taup_model.get_travel_times(
                source_depth_in_km=event['depth'], 
                distance_in_degree=kilometer2degrees(dist/1000.0), 
                phase_list=winspec['phase'].split(','))

        if arrivals:
            f.write('%s|%s|%s|%s,%s,%s|%.1f,%.1f,%.1f|%.1f,%.1f,%.1f|'\
                    '%.4f|%.4f|%.1f|%.1f|%s|%s|%.1f|%.1f|%s|%s\n' % (
                    net_sta_loc[0], net_sta_loc[1], net_sta_loc[2],
                    Z_comp[0][0], H_comp[0][0], H_comp[1][0],
                    Z_comp[0][1], H_comp[0][1], H_comp[1][1],
                    Z_comp[0][2], H_comp[0][2], H_comp[1][2],
                    channels[0]['latitude'], channels[0]['longitude'],
                    channels[0]['elevation'], channels[0]['depth'],
                    winspec['phase'], winspec['orientation'], 
                    winspec['azimuth'], winspec['dip'],
                    event['centroid_time']+arrivals[0].time+winspec['begin'],
                    event['centroid_time']+arrivals[0].time+winspec['end']))
        else:
            print "[WARN] phase not found: ", winspec['phase'], net_sta_loc

# 
f.close()

# END
