from obspy import read_events
from obspy.clients.fdsn import RoutingClient

CMT_file = 'CMTSOLUTION'
evt = read_events(CMT_file)
event = evt[0].origins[0]
origin_time = event.time
starttime = origin_time - 10*60
endtime = origin_time + 30*60
print(event)

client = RoutingClient("iris-federator")

with open('no_station.lst', 'r') as f:
    ns_list = [ l.strip().split('.') for l in f.readlines()]

bulk = [(ns[0], ns[1], '*', '?H?', starttime, endtime) for ns in ns_list]

print(bulk)

inv = client.get_stations_bulk(bulk, latitude=54, longitude=4, maxradius=40, level='response')

# inv = client.get_stations(
#         level="channel",
#         network='1N',
#         station='HAR2',
#         channel='HH?,BH?',
#         starttime=origin_time - 10*60,
#         endtime=origin_time + 30*60,
#         latitude=54,
#         longitude=4,
#         maxradius=40
#         )

inv.write('channel.xml', format='stationxml')

print(inv)
