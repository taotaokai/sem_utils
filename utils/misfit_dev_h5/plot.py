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

args = parser.parse_args()

def plot(arg):
    print(arg)
    misfit = Misfit(args.misfit_h5file)
    misfit.plot_seismogram_1comp(arg[0], out_fig=arg[1])

if __name__ == '__main__':
    print(f"======= [{datetime.now()}] plot_seismogram_1comp")

    with pt.open_file(args.misfit_h5file, 'r') as h5f:
        win_tbl = h5f.root['window']
        win_ids = [w.decode() for w in set(win_tbl.read(field='id'))]
    inputs = [(win_id, os.path.join(args.out_fig_dir, f"{win_id}.pdf")) for win_id in win_ids]

    pool = multiprocessing.Pool(processes=args.nproc)
    pool.map(plot, inputs)

    print(f"======= [{datetime.now()}] END")