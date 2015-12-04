#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import json
import numpy as np

#------
# read command line args
misfit_file = str(sys.argv[1])

#------
#print "\ninitialize\n"
misfit = Misfit()

#------
#print "\nload data\n"
misfit.load(filename=misfit_file)

#------
#print "\nprint misfit\n"

events = misfit.data['events']
for event_id, event in events.iteritems():
    #    print "# %s" % event_id
    stations = event['stations']
    reloc = event['relocate']
    dm = reloc['dm']
    print "#RELOC %-16s dNorth %12.5e dEast %12.5e dDepth %12.5e dT %12.5e" % (
            event_id, dm['dNorth'], dm['dEast'], dm['dDepth'], dm['dT'])
    for station_id, station in stations.iteritems():
        if station['stat']['code'] < 0: 
            #print "#%s %s NOT USED" % (event_id, station_id)
            continue
        windows = station['windows']
        for window_id, window in windows.iteritems():
            if window['stat']['code'] < 0: 
                #print "#%s %s %s NOT USED" % (event_id, station_id, window_id)
                continue
            cc = window['misfit']
            snr = window['quality']['SNR']
            ttime = window['phase']['ttime']
            print "%-16s %-12s %-8s %7.4f %7.4f %8.4f %10.3e %7.3f %7.3f %5.1f %5.3f" % (
                    event_id, station_id, window_id,
                    cc['CC0'],cc['CCmax'],cc['CC_time_shift'],
                    cc['CC_time_shift']/ttime,
                    cc['AR0'],cc['ARmax'],snr, window['weight'])