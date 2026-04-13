#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot misfit statistics from CSV files.

This script reads misfit data from multiple CSV files, processes the data,
and generates a plot showing the weighted cross-correlation sum over iterations.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def adjust_min_gap(arr, min_gap):
    """
    Adjust array values to maintain minimum gap between consecutive elements.
    
    This function sorts the array, ensures each element is at least min_gap
    larger than the previous one, then restores the original order.
    
    Args:
        arr (numpy.ndarray): Input array to adjust
        min_gap (float): Minimum gap between consecutive elements
        
    Returns:
        numpy.ndarray: Adjusted array with maintained minimum gaps
    """
    idx_sort = np.argsort(arr)
    sorted_arr = arr[idx_sort]
    adjusted = np.empty_like(arr)
    adjusted[0] = sorted_arr[0]
    
    for i in range(1, len(arr)):
        new_val = max(sorted_arr[i], adjusted[i-1] + min_gap)
        adjusted[i] = new_val
    
    idx_reverse = np.argsort(idx_sort)
    return adjusted[idx_reverse]

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--fig_title", required=True, help="figure title, e.g. event name")
parser.add_argument("--csv_list", required=True, help="File containing list of CSV files to process")
parser.add_argument("--out_fig", required=True, help="Output figure")
parser.add_argument("--out_csv", default=None, help="Output csv file")
parser.add_argument("--fig_height", default=6, type=float, help="figure height")
parser.add_argument("--fig_width", default=8, type=float, help="figure width")
args = parser.parse_args()

# Read list of CSV files
with open(args.csv_list, "r") as f:
    csv_files = [line.strip() for line in f.readlines()]

# Load and concatenate all CSV data
df = pd.concat([pd.read_csv(csv_file) for csv_file in csv_files], ignore_index=True)

# Create index column combining stage and iteration information
df['index'] = df.apply(lambda row: f"s{row['stage']:02d}_i{row['iter']:02d}", axis=1)
df = df.sort_values(['stage', 'iter'])

# Create mapping from index strings to x-axis positions
indices = df['index'].unique()
xticks = {v: i for i, v in enumerate(indices)}
max_x = len(indices) - 1

# Create the plot
fig = plt.figure(figsize=[args.fig_width, args.fig_height])
ax = fig.add_axes([0.1, 0.1, 0.6, 0.8]) 
ax_height = 0.8 * 8 # in inches
legend_data = []

# Group by window and index, sum over the same (window, index)
df = df.groupby(['window', 'index'])[['wcc0_sum', 'weight_sum']].sum()
df['wcc0_avg'] = df['wcc0_sum'] / df['weight_sum']
if args.out_csv:
    df.to_csv(args.out_csv)

# Plot data for each window
# for win, group in df.groupby('window'):
for win in df.index.get_level_values('window').unique():
    df_sel = df.loc[win]
    # y = group["wcc0_sum"].tolist()
    # x = [xticks[v] for v in group['index']]
    # wcc = group.groupby("index")['wcc0_sum'].sum()  # sum over same index
    # df_sel = df.sort_index()
    y = df_sel['wcc0_sum'].tolist()
    x = [xticks[v] for v in df_sel.index]
    line, = ax.plot(x, y, "-s")
    legend_data.append((y[-1], win, line.get_color()))

# Adjust legend positions to avoid overlap
y_values = np.array([v[0] for v in legend_data])
min_y, max_y = y_values.min(), y_values.max()
text_fontsize = 8
# Calculate minimum gap based on plot height and font size
min_gap = (max_y - min_y) / ax_height * (text_fontsize / 72) * 1.5
adjusted_y = adjust_min_gap(y_values, min_gap)

# Draw connector lines and add legend labels
for i, (original_y, win, color) in enumerate(legend_data):
    ax.plot([max_x, max_x+0.1], [original_y, adjusted_y[i]], clip_on=False, color=color, linewidth=1.0)
    ax.text(max_x+0.1, adjusted_y[i], f" {win}", color=color, fontsize=text_fontsize, va="center", clip_on=False)

# Configure plot appearance
ax.set_ylabel("sum(weight*cc0)")
ax.set_xticks(range(len(indices)))
ax.set_xticklabels(indices)
ax.set_xlim([0, max_x])
ax.grid()
ax.set_title(args.fig_title)
plt.savefig(args.out_fig)