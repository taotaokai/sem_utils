#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" get model update step length from grid search results
"""
import argparse
import tables as pt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )
parser.add_argument("out_file")
parser.add_argument("--phase", default=None)
parser.add_argument("--freqmax", default=None, type=float)
parser.add_argument("--cmpnm", nargs="+", default=['Z', 'R', 'T'])
parser.add_argument("--min_snr", default=5.0, type=float)

args = parser.parse_args()

misfit_h5file = args.misfit_h5file

stats = []

# create table
with pt.open_file(misfit_h5file, "r") as h5f:
    if "/window" not in h5f:
        msg = f"No /window in {misfit_h5file}"
        raise Exception(msg)
    tbl = h5f.get_node("/window")

    tbl_sta = h5f.get_node("/station")

    net_sta_list = set((net, sta) for net, sta in zip(tbl.cols.network, tbl.cols.station))
    net_sta_list = sorted(net_sta_list)
    for net, sta in net_sta_list:
        print(f"[INFO] {net=}, {sta=}")
        stations = tbl_sta.read_where(f"(network == {net}) & (station == {sta})")
        if len(stations) != 1:
            msg = f"{net=}, {sta=} has {len(stations)} stations"
            raise Exception(msg)
        stainfo = stations[0]
        evla = stainfo["latitude"]
        evlo = stainfo["longitude"]
        print(f"[INFO] {net=}, {sta=}, {evla=}, {evlo=}")

        windows = tbl.read_where(f"(network == {net}) & (station == {sta})")
        if args.phase is not None:
            windows = windows[windows['phase'] == args.phase.encode()]
        windows = windows[windows['SNR'] >= args.min_snr]
        nwin = len(windows)
        if nwin == 0:
            continue
        cc0 = windows['cc0']
        ccdt = windows['cc_time_shift']
        snr = windows['SNR']
        win_id = windows['id']
        freq_max = np.array([max(x) for x in windows['butter_Wn']])

        # weight = np.copy(snr)
        # weight[snr > 10] = 10
        # weight[snr < 0] = 0

        data = {}

        data['network'] = net.decode()
        data['station'] = sta.decode()
        data['latitude'] = evla
        data['longitude'] = evlo
        
        for cmpnm in args.cmpnm:
            ind = (windows['cmpnm'] == cmpnm.encode())
            if args.freqmax is not None:
                ind &= (freq_max <= args.freqmax)
            nwin = np.sum(ind)
            data[f"nwin_{cmpnm}"] = nwin
            if nwin == 0:
                continue
            data[f"mean_snr_{cmpnm}"] = np.mean(snr[ind])
            data[f"mean_cc0_{cmpnm}"] = np.mean(cc0[ind]) # * weight[ind]) / np.sum(weight[ind])
            data[f"std_cc0_{cmpnm}"] = np.std(cc0[ind]) # * weight[ind]) / np.sum(weight[ind])
            data[f"mean_ccdt_{cmpnm}"] = np.mean(ccdt[ind])  #* weight[ind]) / np.sum(weight[ind])
            data[f"std_ccdt_{cmpnm}"] = np.std(ccdt[ind])  #* weight[ind]) / np.sum(weight[ind])

        # print(net, sta, win_id, nwin, np.mean(cc0))
        stats.append(data)
    
df = pd.DataFrame(stats)
# df.to_hdf("test.h5", "misfit", mode="w")
df.to_csv(args.out_file, float_format="%.4f", index=False)