import sys
from datetime import datetime
from misfit import Misfit

misfit_h5file   = sys.argv[1]          # "mesh_hdf5/misfit.h5"
out_fig_dir     = sys.argv[2]          # "mesh_hdf5/fig"

with Misfit(misfit_h5file, 'r') as misfit:
    print(f"======= [{datetime.now()}] plot_seismogram_1comp")
    win_tbl = misfit.h5f.root['window']
    win_ids = [w.decode() for w in set(win_tbl.read(field='id'))]
    for win_id in win_ids:
        print(f"[{datetime.now()}] {win_id}")
        misfit.plot_seismogram_1comp(savefig=True, out_dir=out_fig_dir, win_id=win_id)

    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_T_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='p,P_Z_10-50sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_T_10-50sec')
    # # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_Z_30-100sec')
    # # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_R_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='Surf_Z_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='Surf_R_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='Surf_T_30-100sec')

    print(f"======= [{datetime.now()}] END")
