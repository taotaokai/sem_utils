from misfit import Misfit

misfit = Misfit('misfit_test.h5')

misfit.read_config_file('misfit.yaml')

misfit.read_solver_parfile('Par_file')

misfit.read_cmtsolution('CMTSOLUTION.ecef', ECEF=True)

# misfit.write_cmtsolution('C202006211907A.ecef')

misfit.read_channel_file('channel.txt', "STATIONS_FILTERED")

misfit.read_data_h5('sem.h5', "STATIONS_FILTERED")

misfit.read_syn_sac('sac')

misfit.setup_windows()

misfit.measure_adj()

# misfit.output_adj('adj/')

# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='p,P_Z_30-100sec')
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='p,P_Z_10-50sec', align_time=False)
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='s,S_Z_30-100sec')
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='s,S_R_30-100sec')
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='s,S_T_30-100sec')
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='Surf_Z_30-100sec')
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='Surf_R_30-100sec')
# misfit.plot_seismogram_1comp(savefig=True, out_dir='fig/', win_id='Surf_T_30-100sec')

# misfit = Misfit('misfit.h5', 'w')
# misfit.read_config_file('misfit.yaml')
# misfit.close()


# misfit.read_data_h5('data.h5')
# misfit.close()
