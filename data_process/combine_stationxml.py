import argparse

from obspy import read_inventory
from obspy import Inventory

from inventory_utils import isolate_and_merge_station

parser = argparse.ArgumentParser()

parser.add_argument("sxml_list")
parser.add_argument("network")
parser.add_argument("station")
parser.add_argument("out_file")
parser.add_argument("--ignore_station_creation_date", action="store_true")

args = parser.parse_args()
print(args)

# read all xml files
inv_all = Inventory()
with open(args.sxml_list, "r") as f:
    xml_files = [l.strip() for l in f.readlines()]
    for xml_file in xml_files:
        inv = read_inventory(xml_file)
        inv_all.extend(inv.select(network=args.network, station=args.station))

# merge stations and channels
inv = isolate_and_merge_station(
    inv_all,
    args.network,
    args.station,
    ignore_station_creation_date=args.ignore_station_creation_date,
)

# write out station xml
inv.write(args.out_file, format="stationxml")
