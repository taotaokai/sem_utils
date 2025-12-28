import os
from obspy import read, read_inventory, UTCDateTime
from obspy.clients.filesystem.tsindex import Client

waveforms_path = 'waveforms/'

t1 = UTCDateTime('20200321T0125Z')

client = Client('waveforms/timeseries.sqlite', datapath_replace=("^", os.path.join(waveforms_path,'')))

st = client.get_waveforms('NO', 'BRBA', '', 'HHZ',
                          # UTCDateTime('2020-03-21T00:34:53.000000Z'),
                          # UTCDateTime('2020-03-21T01:24:57.463000Z'),
                          # UTCDateTime('2020-03-21T00:35:10.150000Z'),
                          UTCDateTime('2020-03-21T00:35:14.170000Z'),
                          UTCDateTime('2020-03-21T00:35:14.190000Z'),
                          merge=True,
                          )

print(st)
