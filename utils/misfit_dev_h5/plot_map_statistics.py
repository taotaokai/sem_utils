#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" get model update step length from grid search results
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import yaml

parser = argparse.ArgumentParser()

parser.add_argument("misfit_csv" )
parser.add_argument("misfit_par" )
parser.add_argument("out_fig" )
parser.add_argument("--highlight_low_cc0", default=0.0, type=float)

args = parser.parse_args()

column_names = ["mean_cc0_Z", "mean_cc0_R", "mean_cc0_T"]
ncol = len(column_names)

df = pd.read_csv(args.misfit_csv)

df = df.sort_values(column_names, na_position='first')
stlo = df["longitude"].values
stla = df["latitude"].values
net = df["network"].values
sta = df["station"].values

with open(args.misfit_par, 'r') as file:
    config = yaml.safe_load(file)

CRS = getattr(ccrs, config["plot"]["map_type"])
map_params = config["plot"]["map_params"]
projection = CRS(**map_params)
map_extent = config["plot"]["map_extent"]

figsize = (5, 5 * ncol / 2)
fig = plt.figure(figsize=figsize)

gs = GridSpec(ncol+1, 1, height_ratios=[1,] * ncol + [0.05, ],)
              # left=0.05, right=0.95, bottom=0.1, top=0.95,
              #wspace=0.1, hspace=0.1)

axs = []
for i in range(len(column_names)):
    ax = fig.add_subplot(gs[i, 0], projection=projection)
    axs.append(ax)
cax = fig.add_subplot(gs[3, 0])

norm = mcolors.Normalize(vmin=-1, vmax=1)
images = []
for ax, col in zip(axs, column_names):
    ax.set_extent(map_extent)
    ax.coastlines(resolution='50m', color='black', linewidth=0.2, alpha=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.2, alpha=0.5)
    ax.gridlines(draw_labels={"bottom": "x", "left": "y"},
                 formatter_kwargs={'number_format': '.0f',
                                   'degree_symbol': '',
                                   'direction_label': False},
                 xlabel_style={'fontsize': 7},
                 ylabel_style={'fontsize': 7},
                 xpadding=1.5,
                 ypadding=1.5,
                 linewidth=0.2,)

    val = df[col]
    inds = val > args.highlight_low_cc0
    ax.scatter(stlo[inds], stla[inds], c=val[inds], s=4, marker='^', edgecolors=None, linewidth=0,
               cmap='viridis', norm=norm,
               transform=ccrs.Geodetic(), label=col, alpha=0.5)
    inds = val <= args.highlight_low_cc0
    im = ax.scatter(stlo[inds], stla[inds], c=val[inds], s=4, marker='s',
               cmap='viridis', norm=norm,
               transform=ccrs.Geodetic())

    for i, r in df[inds].iterrows():
        lon, lat, net, sta = r['longitude'], r['latitude'], r['network'], r['station']
        ax.text(lon, lat, f" {net}.{sta}", transform=ccrs.Geodetic(),
                horizontalalignment='left', verticalalignment='center',
                fontsize=2)

    ax.legend(loc="upper right")
    images.append(im)

fig.colorbar(images[0], cax=cax, orientation='horizontal', fraction=0.05, label='mean cc0', shrink=0.2)

fig.savefig(args.out_fig)