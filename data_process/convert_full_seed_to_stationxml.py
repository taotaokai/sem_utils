#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import warnings
from obspy import read_inventory, Inventory

parser = argparse.ArgumentParser()

parser.add_argument("seed_list")
parser.add_argument("out_dir")

args = parser.parse_args()
print(args)

# read metadata from all seed files
with open(args.seed_list, "r") as f:
    seed_files = [l.strip() for l in f.readlines()]

inv_all = Inventory()
for seed_file in seed_files:
    try:
        inv = read_inventory(seed_file)
    except Exception as e:
        warnings.warn(f"failed to read metadata from {seed_file}\nError message:\n{e}")
    inv_all.extend(inv)

# write out stationxml files
for network in inv_all:
    net_code = network.code
    for station in network:
        sta_code = station.code
        inv = inv_all.select(network=net_code, station=sta_code)
        inv.write(os.path.join(args.out_dir, f"{net_code}.{sta_code}.xml"), format="stationxml")
