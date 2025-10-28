#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import warnings
from obspy import read, Stream

parser = argparse.ArgumentParser()

parser.add_argument("seed_list")
parser.add_argument("out_dir")

args = parser.parse_args()
print(args)

# read all seed files
with open(args.seed_list, "r") as f:
    seed_files = [l.strip() for l in f.readlines()]

st_all = Stream()
for seed_file in seed_files:
    try:
        st = read(seed_file)
    except Exception as e:
        warnings.warn(f"failed to read {seed_file}\nError message:\n{e}")
    st_all += st

# get list of unique (net, sta)
net_sta_list = set([(tr.stats.network, tr.stats.station) for tr in st_all])

# write out mseed files
for net, sta in net_sta_list:
    st = st_all.select(network=net, station=sta)
    st.write(os.path.join(args.out_dir, f"{net}.{sta}.mseed"), format="mseed")
