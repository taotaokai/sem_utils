import sys, os
import multiprocessing
import argparse
from datetime import datetime
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("out_fig_dir")     # "fig"
parser.add_argument("-n", "--nproc", default=2, type=int)

args = parser.parse_args()

def plot(arg):
    with Misfit(args.misfit_h5file, 'r') as misfit:
        print(arg)
        misfit.plot_seismogram_1comp(arg[0], out_fig=arg[1])

if __name__ == '__main__':
    print(f"======= [{datetime.now()}] plot_seismogram_1comp")

    with Misfit(args.misfit_h5file, 'r') as misfit:
        win_tbl = misfit.h5f.root['window']
        win_ids = [w.decode() for w in set(win_tbl.read(field='id'))]

    inputs = [(win_id, os.path.join(args.out_fig_dir, f"{win_id}.pdf")) for win_id in win_ids]

    pool = multiprocessing.Pool(processes=args.nproc)
    pool.map(plot, inputs)

    print(f"======= [{datetime.now()}] END")

    # with Misfit(args.misfit_h5file, 'r') as misfit:
    #     print(f"======= [{datetime.now()}] plot_seismogram_1comp")
    #     win_tbl = misfit.h5f.root['window']
    #     win_ids = [w.decode() for w in set(win_tbl.read(field='id'))]
    #     for win_id in win_ids:
    #         print(f"[{datetime.now()}] {win_id}")
    #         fig_fn = f"{win_id}.pdf"
    #         misfit.plot_seismogram_1comp(win_id, out_fig=os.path.join(args.out_fig_dir, fig_fn))
    #     print(f"======= [{datetime.now()}] END")
