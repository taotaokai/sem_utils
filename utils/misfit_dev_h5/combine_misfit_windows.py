#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" get model update step length from grid search results
"""
import sys
import argparse
# from datetime import datetime
import tables as pt
# from scipy.interpolate import RegularGridInterpolator
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr

parser = argparse.ArgumentParser()
parser.add_argument("misfit_h5_list" )  # list of misfit.h5 of each events
parser.add_argument("out_h5")
args = parser.parse_args()

with open(args.misfit_h5_list, 'r') as f:
    misfit_h5_list = [ l.strip() for l in f.readlines() ]

# create table
misfit_h5file = misfit_h5_list[0]
with pt.open_file(misfit_h5file, "r") as h5f:
    if "/window" not in h5f:
        msg = f"No /window in {misfit_h5file}"
        raise Exception(msg)
    tbl = h5f.get_node("/window")
    tbl_dtypes = tbl.coldtypes
# add event name column
tbl_dtypes['event_name'] = "S20"
tbl_dtypes = [ (k, v) for k, v in tbl_dtypes.items() ]
out_h5f = pt.open_file(args.out_h5, "w")
tbl_out = out_h5f.create_table("/", "window", np.dtype(tbl_dtypes), "misfit measurements")

# get misfit window from each misfit.h5
for misfit_h5file in misfit_h5_list:
    print(f"[DEBUG]: Working on {misfit_h5file}", file=sys.stderr, flush=True)
    with pt.open_file(misfit_h5file, "r") as h5f:

        if "/window" not in h5f:
            msg = f"No /window in {misfit_h5file}"
            raise Exception(msg)
        tbl_win = h5f.get_node("/window")

        # event info
        if "/source" not in h5f:
            msg = f"No /source in {misfit_h5file}"
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        evnm = event["id"].decode()

        row_list = []
        for row in tbl_win:
            data = tuple(list(row[:]) + [evnm, ])
            row_list.append(data)

        tbl_out.append(row_list)
        tbl_out.flush()

out_h5f.close()