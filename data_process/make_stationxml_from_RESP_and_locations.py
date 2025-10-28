# import sys
import os
import warnings
# import copy
import argparse

import pandas as pd
from obspy import read_inventory
# from obspy.core.inventory.util import Equipment
from obspy import Inventory

from inventory_utils import isolate_and_merge_station

parser = argparse.ArgumentParser()

# def _to_datetime(value):
#     try:
#         datetime = pd.Timestamp(value)
#     except Exception as e:
#         raise argparse.ArgumentTypeError(f"{value} is not a valid datetime string! ({e})")
#     return datetime

parser.add_argument("resp_list")
# parser.add_argument("event_date", type=_to_datetime)
parser.add_argument("station_file")
parser.add_argument("out_dir")

args = parser.parse_args()
print(args)

# event_date = args.event_date
station_file = args.station_file
out_dir = args.out_dir

# # metadata_df = pd.read_excel(metadata_xls)
# metadata_df = pd.read_csv(metadata_csv)
# # filter metadata by event_date
# metadata_df['In Date'] = pd.to_datetime(metadata_df['In Date'], format="%Y%m%d")
# metadata_df.loc[metadata_df['Out Date'] == "OnSite", 'Out Date'] = "20990101"
# metadata_df['Out Date'] = pd.to_datetime(metadata_df['Out Date'], format="%Y%m%d")
# mask = (metadata_df['In Date'] <= event_date) & (metadata_df['Out Date'] >= event_date)
# metadata_df = metadata_df[mask]

# read station file
metadata_df = pd.read_csv(station_file, sep=r"\s+", header=None)
metadata_df.columns = ['net', 'sta', 'lat', 'lon']

# read all RESP files
with open(args.resp_list, "r") as f:
    resp_files = [l.strip() for l in f.readlines()]

inv_all = Inventory()
for resp_file in resp_files:
    try:
        inv = read_inventory(resp_file)
    except Exception as e:
        warnings.warn(f"failed to read invetory file {resp_file}\nError message:\n{e}")
    inv_all.extend(inv)

# get list of net,sta
net_sta_set = set()
for net in inv_all:
    net_code = net.code
    for sta in net:
        sta_code = sta.code
        net_sta_set.add((net_code, sta_code))

# add station metadata
for net_code, sta_code in net_sta_set:
    # print("============ ", net_code, sta_code)
    inv = isolate_and_merge_station(inv_all, net_code, sta_code) # inv_all.select(network=net_code, station=sta_code)

    # get station metadata
    mask = (metadata_df['net'] == net_code) & (metadata_df['sta'] == sta_code)
    sta_df = metadata_df[mask]
    # print(sta_df)
    if len(sta_df) != 1:
        warnings.warn(f"none or more than one station metadata found for {net_code}.{sta_code}")
        # print(cha_df)
        continue
    info = sta_df.iloc[0]
    stla, stlo = info['lat'], info['lon']

    # add station metadata to each channel
    for net in inv:
        for sta in net:
            sta.latitude = stla
            sta.longitude = stlo
            sta.elevation = 0
            channels = []
            for cha in sta:
                cha_code = cha.code
                # get channel metadata
                # mask = sta_df['CMP'] == cha_code
                # cha_df = sta_df[mask]
                # if len(cha_df) != 1:
                #     warnings.warn(f"none or more than one channel metadata found: {net_code}.{sta_code}.{cha_code}")
                #     # print(cha_df)
                #     continue
                # info = cha_df.iloc[0]
                # set channel metadata
                # cha.latitude = info['Latitude (deg)']
                # cha.longitude = info['Longitude (deg)']
                # cha.elevation = info['Elv (m)']
                # description = "%s + %s"%(info['Sensor'], info['Digitizer'])
                # cha.sensor = Equipment(description=description)
                cha.latitude = stla
                cha.longitude = stlo
                cha.elevation = 0
                cha.depth = 0
                # set other channel metadata
                cha.depth = 0.0
                if cha_code[-1] == "E":
                    cha.dip = 0
                    cha.azimuth = 90
                elif cha_code[-1] == "N":
                    cha.dip = 0
                    cha.azimuth = 0
                elif cha_code[-1] == "Z":
                    cha.dip = -90
                    cha.azimuth = 0
                else:
                    warnings.warn(f"unknown channel code: {cha_code}")
                    continue
                channels.append(cha)
            sta.channels = channels

    inv.write(os.path.join(out_dir, f"{net_code}.{sta_code}.xml"), format="stationxml")