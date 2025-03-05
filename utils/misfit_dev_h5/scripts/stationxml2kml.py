import sys
# import glob
from obspy import read_inventory
import simplekml

sxml_file = sys.argv[1]
kml_file = sys.argv[2]
kml_name = sys.argv[3]

# xml_list = glob.glob('stations/*.xml')
# inv_all = Inventory()
# for xml in xml_list:
#     try:
#         inv = read_inventory(xml)
#         inv_all.extend(inv)
#     except Exception as err:
#         print(f'cannot read {xml}!')
#         print(err)
# inv_all.write('stations1.xml', format='stationxml')

inv = read_inventory(sxml_file)

inv_Zcomp = inv.select(channel='[BH]HZ')
channels_Zcomp = inv_Zcomp.get_contents()['channels']

kml = simplekml.Kml()
kml.document.name = kml_name
for chan in channels_Zcomp:
    coord = inv_Zcomp.get_coordinates(chan)
    # coord_info[chan] = coord
    pnt = kml.newpoint(description=chan, coords=[(coord['longitude'], coord['latitude'])])  # lon, lat, optional height

kml.save(kml_file)
