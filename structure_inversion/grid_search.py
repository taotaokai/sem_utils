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
parser.add_argument("--out_nc", default="grid_search_structure.nc")
parser.add_argument("--out_figure", default="grid_search_structure.pdf")
parser.add_argument("--out_txt", default="grid_search_structure.txt")
# parser.add_argument("--step_size", default=0.1, type=float, help="interpolation step size of the grid search results")

args = parser.parse_args()

# sum up grid search results of every events
with open(args.misfit_h5_list, 'r') as f:
    misfit_h5_list = [ l.strip() for l in f.readlines() ]

dm = None
wcc_sum = None
weight_sum = 0.0

for misfit_h5file in misfit_h5_list:
    print(f"[DEBUG]: Working on {misfit_h5file}", file=sys.stderr, flush=True)
    with pt.open_file(misfit_h5file, "r") as h5f:

        if "/grid_search/structure/wcc_sum" not in h5f:
            msg = f"No /grid_search/structure/wcc_sum not in {misfit_h5file}"
            raise Exception(msg)
        data = h5f.get_node("/grid_search/structure/wcc_sum")

        dm_tags = data.attrs['dims']
        weight = data.attrs['weight_sum']
        weight_sum += weight

        wcc = data[:]
        if wcc_sum is None:
            wcc_sum = wcc
        else:
            if wcc_sum.shape != wcc.shape:
                msg = f"{wcc.shape=} in {misfit_h5file} should be equal to {wcc_sum.shape=}"
                raise Exception(msg)
            wcc_sum += wcc

        dm1 = {}
        for tag in dm_tags:
            if f"/grid_search/structure/{tag}" not in h5f:
                msg = f"No /grid_search/structure/{tag} not in {misfit_h5file}"
                raise Exception(msg)
            data = h5f.get_node(f"/grid_search/structure/{tag}")
            dm1[tag] = data[:]
        if dm is None:
            dm = dm1
        else:
            if set(dm.keys()) != set(dm1.keys()) or not all(np.array_equal(dm[k], dm1[k]) for k in dm):
                msg = f"{dm1=} in {misfit_h5file=} should be equal to {dm=}"
                raise Exception(msg)


# # interpolate grid search results onto finer grid
# interp = RegularGridInterpolator(tuple(dm.values()), wcc_sum, method='cubic')
# refined_dm = {tag: np.linspace(dm[tag].min(), dm[tag].max(), n_refined) for tag in dm}


# find the optimal step lengths along each dm
opt_ind = np.unravel_index(np.argmax(wcc_sum, axis=None), wcc_sum.shape)
opt_dm = [dm[tag][opt_ind[i]] for i, tag in enumerate(dm)] 

# save results 
da = xr.DataArray(wcc_sum, dims=list(dm.keys()), coords=dm, name="wcc_sum")
da.attrs['weight_sum'] = weight_sum
da.attrs['opt_ind'] = opt_ind
da.attrs['opt_dm'] = opt_dm
da.to_netcdf(args.out_nc)

# output values of optimal dm
with open(args.out_txt, "w") as f:
    for i, tag in enumerate(dm):
        f.write(f"{tag}  {opt_dm[i]}\n")


# plot profiles of wcc along each dm
n_dm = len(dm)
fig, axs = plt.subplots(n_dm, 1, figsize=(10, 4*n_dm+1))
for i, tag in enumerate(dm):
    # indices along i-th dimension of wcc_sum through optimal index opt_ind
    ind = list(opt_ind)
    ind[i] = slice(None) # slice all values along i-th dimension
    ind = tuple(ind)
    x, y = dm[tag], wcc_sum[ind] # dm and wcc along i-th dimension of wcc_sum
    io = opt_ind[i] # get optimal index on this profile
    xo, yo = x[io], y[io]

    if n_dm == 1:
        ax = axs
    else:
        ax = axs[i]
    ax.plot(x, y, 'ko-')
    ax.plot(xo, yo, 'rs')
    ax.text(xo, yo, f"{tag} = {xo:.2f}")
    ax.set_xlabel(tag)
    ax.set_ylabel("sum(wcc)")
fig.tight_layout()
fig.savefig(args.out_figure) 