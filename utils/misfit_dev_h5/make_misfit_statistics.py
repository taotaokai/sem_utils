#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import tables as pt
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("misfit_h5file" )
parser.add_argument("out_csv")
parser.add_argument("--type", default="source")
parser.add_argument("--stage", default=0, type=int)
parser.add_argument("--iter", default=0, type=int)
args = parser.parse_args()

with pt.open_file(args.misfit_h5file, "r") as h5f:
    if "/window" not in h5f:
        msg = f"No /window in {args.misfit_h5file}"
        raise Exception(msg)
    tbl = h5f.get_node("/window")

    if "/source" not in h5f:
        msg = '"/source" not existing, run read_cmtsolution first!'
        raise KeyError(msg)
    tbl_src = h5f.get_node("/source")
    if tbl_src.nrows == 0:
        msg = "no source information"
        raise Exception(msg)
    event = tbl_src[0]
    evnm = event["id"].decode()

    stats = []
    for win_id in sorted(set(tbl.cols.id)):
        windows = tbl.read_where(f"id == {win_id}")
        cc0 = windows['cc0']
        weights = windows['weight']

        weight_sum = sum(weights)
        wcc0_sum = sum(weights * cc0)
        wcc0_avg = wcc0_sum / weight_sum

        data = {}
        data['type'] = args.type
        data['stage'] = args.stage
        data['iter'] = args.iter
        data['window'] = win_id.decode()
        data["weight_sum"] = weight_sum
        data["wcc0_sum"] = wcc0_sum
        data["wcc0_avg"] = wcc0_avg
        data["event"] = evnm

        # print(net, sta, win_id, nwin, np.mean(cc0))
        stats.append(data)
    
    df = pd.DataFrame(stats)
    df.to_csv(args.out_csv, float_format="%.4f", index=False)