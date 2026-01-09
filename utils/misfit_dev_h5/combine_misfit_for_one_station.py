#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""get model update step length from grid search results"""
import sys
import argparse
import tables as pt
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5_list")  # list of misfit.h5 of each events
parser.add_argument("network", help="network name")
parser.add_argument("station", help="station name")
parser.add_argument("out_h5", help="output h5 file")

args = parser.parse_args()

# sum up grid search results of every events
with open(args.misfit_h5_list, "r") as f:
    misfit_h5_list = [l.strip() for l in f.readlines()]

# create table
misfit_h5file = misfit_h5_list[0]
with pt.open_file(misfit_h5file, "r") as h5f:
    # if "/source" not in h5f:
    #     msg = f"No /source in {misfit_h5file}"
    #     raise Exception(msg)
    # tbl = h5f.get_node("/source")
    # tbl_src_dtypes = tbl.coldtypes
    # tbl_src_dtypes = [(k, v) for k, v in tbl_src_dtypes.items()]
    if "config" not in h5f.root._v_attrs:
        msg = f"No config in root._v_attrs fo {misfit_h5file}"
        raise Exception(msg)
    config = h5f.root._v_attrs["config"]

    if "/window" not in h5f:
        msg = f"No /window in {misfit_h5file}"
        raise Exception(msg)
    tbl = h5f.get_node("/window")
    tbl_win_dtypes = tbl.coldtypes
    tbl_win_dtypes["event_name"] = "S20"
    tbl_win_dtypes = [(k, v) for k, v in tbl_win_dtypes.items()]

out_h5f = pt.open_file(args.out_h5, "w")
out_h5f.root._v_attrs["config"] = config
# out_tbl_src = out_h5f.create_table("/", "source", np.dtype(tbl_src_dtypes), "events")
out_tbl_win = out_h5f.create_table("/", "window", np.dtype(tbl_win_dtypes), "misfit")

net = args.network
sta = args.station

out_h5f.root._v_attrs["network"] = net
out_h5f.root._v_attrs["station"] = sta

# get misfit window from each misfit.h5
for misfit_h5file in misfit_h5_list:
    print(f"[DEBUG]: Working on {misfit_h5file}", file=sys.stderr, flush=True)
    with pt.open_file(misfit_h5file, "r") as h5f:

        # event info
        if "/source" not in h5f:
            msg = f"No /source in {misfit_h5file}"
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        # out_tbl_src.append([tuple(event), ])
        # out_tbl_src.flush()

        evnm = event["id"].decode()
        evlo = event["longitude"]
        evla = event["latitude"]
        evdp = event["depth"]
        evmt = event["mt_rtp"]
        evt0 = event["t0"]
        evtau = event["tau"]

        if "/window" not in h5f:
            msg = f"No /window in {misfit_h5file}"
            raise Exception(msg)
        tbl_win = h5f.get_node("/window")

        row_list = []
        for row in tbl_win.where(f'(network == b"{net}") & (station == b"{sta}")'):
            data = tuple(list(row[:]) + [evnm])
            row_list.append(data)
        if len(row_list) == 0:
            continue
        out_tbl_win.append(row_list)
        out_tbl_win.flush()

        path_sta = f"/waveform/{net}_{sta}"
        if path_sta not in h5f:
            msg = f"No {path_sta} in {misfit_h5file}"
            raise KeyError(msg)
        grp_sta = h5f.get_node(path_sta)
        stla = grp_sta._v_attrs["latitude"]
        stlo = grp_sta._v_attrs["longitude"]
        stel = grp_sta._v_attrs["elevation"]

        out_grp = out_h5f.create_group("/waveform", f"{evnm}", createparents=True)
        out_grp._v_attrs["event_name"] = evnm
        out_grp._v_attrs["origin_time"] = evt0
        out_grp._v_attrs["evla"] = evla
        out_grp._v_attrs["evlo"] = evlo
        out_grp._v_attrs["evdp"] = evdp
        out_grp._v_attrs["evtau"] = evtau
        out_grp._v_attrs["moment_tensor"] = evmt
        out_grp._v_attrs["stla"] = stla
        out_grp._v_attrs["stlo"] = stlo
        out_grp._v_attrs["stel"] = stel

        # Copy the node to the destination file's root group
        if "DATA_DISP" in grp_sta:
            array_node = grp_sta["DATA_DISP"]
            h5f.copy_node(array_node,
                          newparent=out_grp, 
                          overwrite=True)
        if "DATA_VEL" in grp_sta:
            array_node = grp_sta["DATA_VEL"]
            h5f.copy_node(array_node,
                          newparent=out_grp, 
                          overwrite=True)
        if "SYN_DISP" in grp_sta:
            array_node = grp_sta["SYN_DISP"]
            h5f.copy_node(array_node,
                          newparent=out_grp, 
                          overwrite=True)

out_h5f.close()
