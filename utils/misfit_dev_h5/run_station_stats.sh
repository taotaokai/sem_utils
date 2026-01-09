find events -wholename "*/misfit/misfit.h5" > event_misfit_h5.lst

python sem_utils/misfit/combine_misfit_windows.py event_misfit_h5.lst misfit_windows.h5

python sem_utils/misfit/make_station_statistics.py misfit_windows.h5 misfit_surf_0.05Hz.csv --phase surf --freqmax 0.055 --min_snr 5

python sem_utils/misfit/plot_station_statistics.py misfit_surf_0.05Hz.csv misfit_surf_0.05Hz.pdf
