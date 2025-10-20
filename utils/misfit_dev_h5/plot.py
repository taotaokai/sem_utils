import sys, os
import multiprocessing
import argparse
from datetime import datetime
import tables as pt
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("out_fig_dir")     # "fig"
parser.add_argument("-n", "--nproc", default=2, type=int)
parser.add_argument("--title", default=None)

args = parser.parse_args()

def plot(arg):
    print(arg)
    misfit = Misfit(args.misfit_h5file)
    misfit.plot_seismogram_1comp(arg[0], out_fig=arg[1], title_prefix=args.title)

if __name__ == '__main__':
    print(f"======= [{datetime.now()}] plot_seismogram_1comp")

    misfit = Misfit(args.misfit_h5file)

    pool = multiprocessing.Pool(processes=args.nproc)
    window_ids = misfit.get_window_ids()
    inputs = [(win_id, os.path.join(args.out_fig_dir, f"{win_id}.pdf")) for win_id in window_ids]
    pool.map(plot, inputs)

    print(f"======= [{datetime.now()}] plot_histogram")
    misfit.plot_histogram(out_fig=os.path.join(args.out_fig_dir, "misfit_hist.pdf"), title_prefix=args.title)

    print(f"======= [{datetime.now()}] END")
