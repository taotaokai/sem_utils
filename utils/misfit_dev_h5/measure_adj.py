import sys
from datetime import datetime
from misfit import Misfit

misfit_parfile  = sys.argv[1]          # "misfit.yaml"
solver_parfile  = sys.argv[2]          # "mesh_hdf5/DATA/Par_file"
cmt_file        = sys.argv[3]          # "mesh_hdf5/DATA/CMTSOLUTION"
channel_file    = sys.argv[4]          # "data/C202006211907A/channel.txt"
data_h5file     = sys.argv[5]          # "data/C202006211907A/data.h5"
syn_sac_dir     = sys.argv[6]          # "mesh_hdf5/OUTPUT_FILES/sac"

misfit_h5file   = sys.argv[7]          # "mesh_hdf5/misfit.h5"
out_adj_dir     = sys.argv[8]          # "mesh_hdf5/adj"
syn_is_grn      = sys.argv[9]          # "Ture"
# out_fig_dir       = sys.argv[2]          # "mesh_hdf5/fig"

is_grn = False
if syn_is_grn == "True":
  is_grn = True

with Misfit(misfit_h5file, 'w') as misfit:
    print(f"======= [{datetime.now()}] read_config_file")
    misfit.read_config_file(misfit_parfile)
    print(f"======= [{datetime.now()}] read_solver_parfile")
    misfit.read_solver_parfile(solver_parfile)
    print(f"======= [{datetime.now()}] read_cmtsolution")
    misfit.read_cmtsolution(cmt_file, ECEF=True)
    print(f"======= [{datetime.now()}] read_channel_file")
    misfit.read_channel_file(channel_file)
    print(f"======= [{datetime.now()}] read_data_h5")
    misfit.read_data_h5(data_h5file)
    print(f"======= [{datetime.now()}] read_syn_sac")
    misfit.read_syn_sac(syn_sac_dir, is_grn=is_grn)
    print(f"======= [{datetime.now()}] setup_windows")
    misfit.setup_windows()
    print(f"======= [{datetime.now()}] measure_adj")
    misfit.measure_adj()
    print(f"======= [{datetime.now()}] output_adj")
    misfit.output_adj(out_adj_dir)
    # print(f"======= [{datetime.now()}] plot_seismogram_1comp")
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='p,P_Z_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_T_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='p,P_Z_10-50sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_T_10-50sec')
    # # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_Z_30-100sec')
    # # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='s,S_R_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='Surf_Z_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='Surf_R_30-100sec')
    # misfit.plot_seismogram_1comp(savefig=True, out_dir=fig_dir, win_id='Surf_T_30-100sec')
