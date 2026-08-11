#!/usr/bin/env python
"""
Plot misfit statistics from CSV files.

This script reads misfit data from multiple CSV files, processes the data,
and generates a plot showing the weighted cross-correlation sum over iterations.
"""

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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
        new_val = max(sorted_arr[i], adjusted[i - 1] + min_gap)
        adjusted[i] = new_val

    idx_reverse = np.argsort(idx_sort)
    return adjusted[idx_reverse]


# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--fig_title", required=True, help="figure title, e.g. event name")
parser.add_argument(
    "--csv_list", required=True, help="File containing list of CSV files to process"
)
parser.add_argument("--out_fig", required=True, help="Output figure")
parser.add_argument("--fig_height", default=6, type=float, help="figure height")
parser.add_argument("--fig_width", default=8, type=float, help="figure width")
args = parser.parse_args()

# Read list of CSV files
with open(args.csv_list, "r") as f:
    csv_files = [line.strip() for line in f]

# Load and concatenate all CSV data
df = pd.concat([pd.read_csv(csv_file) for csv_file in csv_files], ignore_index=True)

# Create index column combining stage and iteration information
df['stage_iter'] = df.apply(lambda row: f"s{row['stage']:02d}_i{row['iter']:02d}", axis=1)
df = df.sort_values(['stage', 'iter', 'event'])

# Create mapping from event strings to x-axis positions
indices = df["event"].unique()
xticks = {v: i for i, v in enumerate(indices)}
max_x = len(indices) - 1

# Group by stage_iter,event
df_all = df.groupby(["stage_iter", "event"])[["wcc0_sum", "weight_sum"]].sum()
df_all["wcc0_avg"] = df_all["wcc0_sum"] / df_all["weight_sum"]

# Group by window and event, sum over the same (window, event)
df = df.groupby(["window", "stage_iter", "event"])[["wcc0_sum", "weight_sum"]].sum()
df["wcc0_avg"] = df["wcc0_sum"] / df["weight_sum"]

min_wcc0_avg = df["wcc0_avg"].min()
max_wcc0_avg = df["wcc0_avg"].max()
min_wcc0_sum = df["wcc0_sum"].min()
max_wcc0_sum = df["wcc0_sum"].max()


def plot_data(df, fig_title, min_wcc0_sum, max_wcc0_sum):
    # Create the plot
    fig = plt.figure(figsize=[args.fig_width, args.fig_height])
    ax = fig.add_axes([0.1, 0.15, 0.8, 0.8])
    ax_height = 0.8 * args.fig_height  # in inches
    ax2 = ax.twinx()

    legend_data = []

    stage_iter_list = df.index.get_level_values("stage_iter").unique()
    niter = len(stage_iter_list)
    width = 0.8 / niter
    for i, stage_iter in enumerate(stage_iter_list):
        df_stage = df.loc[stage_iter]
        y = df_stage["wcc0_avg"].tolist()
        x = [xticks[v] for v in df_stage.index]
        (line,) = ax.plot(x, y, "-o")
        line_color = line.get_color()
        legend_data.append((y[-1], stage_iter, line_color))

        y = df_stage["wcc0_sum"].tolist()
        x = np.array([xticks[v] for v in df_stage.index])
        offset = -0.4 + i * width
        ax2.bar(x + offset, y, width=width, align='edge', color=line_color, alpha=0.4)
        # legend_data.append((y[-1], stage_iter, line.get_color()))

    # y = group["wcc0_sum"].tolist()
    # x = [xticks[v] for v in group['index']]
    # wcc = group.groupby("index")['wcc0_sum'].sum()  # sum over same index
    # df_sel = df.sort_index()
    # y = df_sel["wcc0_avg"].tolist()
    # x = [xticks[v] for v in df_sel.index]
    # (line,) = ax.plot(x, y, "-ko")

    # z = df_sel['wcc0_sum'].tolist()
    # for i in range(len(z)):
    #     ax.text(x[i], y[i], f"{z[i]:.1f}", va="top", ha="center")

    # legend_data.append((y[-1], win, line.get_color()))

    # Adjust legend positions to avoid overlap
    y_values = np.array([v[0] for v in legend_data])
    text_fontsize = 18
    # Calculate minimum gap based on plot height and font size
    min_gap = (max_wcc0_avg - min_wcc0_avg) / ax_height * (text_fontsize / 72) * 1.5
    adjusted_y = adjust_min_gap(y_values, min_gap)
    
    # Draw connector lines and add legend labels
    for i, (original_y, stage_iter, color) in enumerate(legend_data):
        ax.plot(
            [max_x, max_x + 1.5],
            [original_y, adjusted_y[i]],
            clip_on=False,
            color=color,
            linewidth=1.0,
        )
        ax.text(
            max_x + 1.5,
            adjusted_y[i],
            f" {stage_iter}",
            color=color,
            fontsize=text_fontsize,
            va="center",
            clip_on=False,
        )

    # Configure plot appearance
    ax.set_ylabel("weighted average of cc0")
    ax.set_ylim([min_wcc0_avg-0.1, max_wcc0_avg+0.1])
    ax.set_xticks(range(len(indices)))
    ax.set_xticklabels(indices, rotation=90, ha="center")
    ax.set_xlim([-1, max_x + 1])
    ax.grid()
    # ax.set_title(f"{args.fig_title}: {win}")
    ax.set_title(fig_title)

    ax2.set_ylabel("sum(weight*cc0)")
    # ax2.set_ylim([df['wcc0_sum'].min(), df['wcc0_sum'].max()])
    ax2.set_ylim([min_wcc0_sum, max_wcc0_sum])

    return fig



# Plot data for each window
# for win, group in df.groupby('window'):
with PdfPages(args.out_fig) as pdf:

    fig_title = f"{args.fig_title}: all"
    fig = plot_data(df_all, fig_title, df_all["wcc0_sum"].min(), df_all["wcc0_sum"].max()) 
    pdf.savefig(fig)
    plt.close(fig)

    for win in df.index.get_level_values("window").unique():
        fig_title = f"{args.fig_title}: {win}"
        df_win = df.loc[win]
        fig = plot_data(df_win, fig_title, min_wcc0_sum, max_wcc0_sum) 
        pdf.savefig(fig)
        plt.close(fig)

    

