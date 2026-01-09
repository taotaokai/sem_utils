#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" get model update step length from grid search results
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("misfit_csv" )
parser.add_argument("out_fig" )

args = parser.parse_args()

df = pd.read_csv(args.misfit_csv)

df = df.sort_values(["mean_cc0_Z", "mean_cc0_R", "mean_cc0_T"], na_position='first')
nsta = df.shape[0]

figsize = (5, nsta//10)
fig, ax = plt.subplots(figsize=figsize)

fig.subplots_adjust(
    left=0.10,
    right=0.85,
    bottom=0.02,
    top=0.98,
    # wspace=0.1,  # Reduced horizontal space
    # hspace=0.1   # Reduced vertical space
)

# ax.tick_params(
#     axis='both',
#     bottom=True,
#     top=True,
#     left=True,
#     right=True,
#     labelbottom=True,
#     labeltop=True,
#     labelleft=True,
#     labelright=True
# )

y = np.arange(nsta)
ax.scatter(df['mean_cc0_Z'], y, s=5, facecolor='k', marker='o', label="Z")
ax.scatter(df['mean_cc0_R'], y, s=5, facecolor='b', marker='o', label="R")
ax.scatter(df['mean_cc0_T'], y, s=5, facecolor='r', marker='^', label="T")

mean, std = df["mean_cc0_Z"], df["std_cc0_Z"]
ax.hlines(y=y,xmin=mean-std,xmax=mean+std, linewidth=5.0, color='k', alpha=0.2)
mean, std = df["mean_cc0_R"], df["std_cc0_R"]
ax.hlines(y=y,xmin=mean-std,xmax=mean+std, linewidth=2.0, color='b', alpha=0.5)
mean, std = df["mean_cc0_T"], df["std_cc0_T"]
ax.hlines(y=y,xmin=mean-std,xmax=mean+std, linewidth=0.5, color='r', alpha=1.0)

for i in range(nsta):
    info = df.iloc[i]

    annot_Z = f"{info["mean_snr_Z"]:.0f}({info["nwin_Z"]})"
    annot_R = f"{info["mean_snr_R"]:.0f}({info["nwin_R"]})"
    annot_T = f"{info["mean_snr_T"]:.0f}({info["nwin_T"]})"
    ax.text(1, y[i], f"{annot_Z},{annot_R},{annot_T}", fontsize=5, zorder=3)

    annot_Z = f"{info["mean_ccdt_Z"]:.1f}({info["std_ccdt_Z"]})"
    annot_R = f"{info["mean_ccdt_R"]:.1f}({info["std_ccdt_R"]})"
    annot_T = f"{info["mean_ccdt_T"]:.1f}({info["std_ccdt_T"]})"
    ax.text(-1, y[i], f"{annot_Z},{annot_R},{annot_T}", fontsize=5, zorder=3)


ax.set_xticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
station_labels = [ f"{net}.{sta}" for net, sta in zip(df['network'].values, df['station'].values) ]
ax.set_yticks(ticks=y, labels=station_labels, fontsize=5)

ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1, nsta)
ax.grid(linewidth=0.2, alpha=0.5)
ax.vlines(x=0, ymin=-1, ymax=nsta, linewidth=0.5, color='k')

ax.set_xlabel("mean CC0")
ax.legend(loc="lower right")

fig.savefig(args.out_fig)