import sys, os
import argparse
from datetime import datetime
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("out_fig_dir")     # "fig"

args = parser.parse_args()

with Misfit(args.misfit_h5file, 'r') as misfit:
    print(f"======= [{datetime.now()}] plot_seismogram_1comp")
    win_tbl = misfit.h5f.root['window']
    win_ids = [w.decode() for w in set(win_tbl.read(field='id'))]
    for win_id in win_ids:
        print(f"[{datetime.now()}] {win_id}")
        fig_fn = f"{win_id}.pdf"
        misfit.plot_seismogram_1comp(win_id, out_fig=os.path.join(args.out_fig_dir, fig_fn))
    print(f"======= [{datetime.now()}] END")
