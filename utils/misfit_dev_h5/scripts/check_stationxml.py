import glob
from obspy import read, read_inventory, Inventory
import simplekml

xml_list = glob.glob('stations/*.xml')

inv_all = Inventory()
for xml in xml_list:
    try:
        inv = read_inventory(xml)
        inv_all.extend(inv)
    except Exception as err:
        print(f'cannot read {xml}!')
        print(err)
inv_all.write('stations1.xml', format='stationxml')

# channels = inv_all.get_contents()['channels']
# channels_Zcomp = [chan for chan in channels if chan[-1]=='Z']
inv_Zcomp = inv_all.select(channel='[BH]HZ')
channels_Zcomp = inv_Zcomp.get_contents()['channels']

kml = simplekml.Kml()
kml.document.name = '202003210049A'
for chan in channels_Zcomp:
    coord = inv_Zcomp.get_coordinates(chan)
    # coord_info[chan] = coord
    pnt = kml.newpoint(description=chan, coords=[(coord['longitude'], coord['latitude'])])  # lon, lat, optional height

kml.save('stations1.kml')
