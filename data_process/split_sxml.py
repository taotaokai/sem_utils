#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import warnings
from obspy import read_inventory, Inventory

parser = argparse.ArgumentParser()

parser.add_argument("sxml_list")
parser.add_argument("out_dir")

args = parser.parse_args()
print(args)

# read metadata from all seed files
with open(args.sxml_list, "r") as f:
    sxml_files = [l.strip() for l in f.readlines()]

for sxml_file in sxml_files:
    try:
        inv = read_inventory(sxml_file)
    except Exception as e:
        warnings.warn(f"failed to read metadata from {sxml_file}\nError message:\n{e}")

    # write out individual stationxml files for each station
    for network in inv:
        net_code = network.code
        for station in network:
            sta_code = station.code
            inv1 = inv.select(network=net_code, station=sta_code)
            out_path = os.path.join(args.out_dir, f"{net_code}.{sta_code}.xml")
            if os.path.exists(out_path):
                msg = f"{out_path=} exists, overwriting!!!"
                warnings.warn(msg)
            inv1.write(out_path, format="stationxml")

