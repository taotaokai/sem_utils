# from obspy.clients.fdsn import RoutingClient
# client = RoutingClient("eida-routing")
# starttime = '2020-03-21T00:39:54.30'
# endtime = '2020-03-21T01:19:54.30'
# inv = client.get_stations(starttime=starttime, endtime=endtime, latitude=54, longitude=4, maxradius=40, channel='?HZ', level='channel')
# print(inv)

from obspy import read_events
from obspy.clients.fdsn.mass_downloader import CircularDomain, \
    Restrictions, MassDownloader

CMT_file = 'CMTSOLUTION'
evt = read_events(CMT_file)
event = evt[0].origins[0]
origin_time = event.time
print(event)

# Circular domain around the epicenter. This will download all data between
# 70 and 90 degrees distance from the epicenter. This module also offers
# rectangular and global domains. More complex domains can be defined by
# inheriting from the Domain class.
domain = CircularDomain(latitude=54, longitude=4, minradius=0.0, maxradius=40.0)

restrictions = Restrictions(
    # Get data from 5 minutes before the event to one hour after the
    # event. This defines the temporal bounds of the waveform data.
    starttime=origin_time - 10*60,
    endtime=origin_time + 30*60,
    # You might not want to deal with gaps in the data. If this setting is
    # True, any trace with a gap/overlap will be discarded.
    reject_channels_with_gaps=True,
    # And you might only want waveforms that have data for at least 95 % of
    # the requested time span. Any trace that is shorter than 95 % of the
    # desired total duration will be discarded.
    minimum_length=0.95,
    # No two stations should be closer than 1 km to each other. This is
    # useful to for example filter out stations that are part of different
    # networks but at the same physical station. Settings this option to
    # zero or None will disable that filtering.
    minimum_interstation_distance_in_m=1E3,
    # Only HH or BH channels. If a station has HH channels, those will be
    # downloaded, otherwise the BH. Nothing will be downloaded if it has
    # neither. You can add more/less patterns if you like.
    channel_priorities=["BH[ZNE12]", "HH[ZNE12]"],
    # Location codes are arbitrary and there is no rule as to which
    # location is best. Same logic as for the previous setting.
    # location_priorities=["", "00", "10"]
    )

# No specified providers will result in all known ones being queried.
mdl = MassDownloader()
# The data will be downloaded to the ``./waveforms/`` and ``./stations/``
# folders with automatically chosen file names.
mdl.download(domain, restrictions, mseed_storage="waveforms",
             stationxml_storage="stations")

