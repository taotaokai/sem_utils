import sys
import pandas as pd
import numpy as np
from obspy import UTCDateTime

sem_parfile = sys.argv[1]

sem_params = pd.read_csv(
    sem_parfile,
    delimiter=r"\s*=\s*",
    header=None,
    comment="#",
    names=["key", "value"],
    dtype=dict(key=object, value=object),
    index_col=["key"],
    engine='python',
).to_dict()["value"]

mesh_central_lat = sem_params['CENTER_LATITUDE_IN_DEGREES']
mesh_central_lon = sem_params['CENTER_LONGITUDE_IN_DEGREES']
mesh_width_xi =    sem_params['ANGULAR_WIDTH_XI_IN_DEGREES']
mesh_width_eta =   sem_params['ANGULAR_WIDTH_ETA_IN_DEGREES']
mesh_gamma_rot =   sem_params['GAMMA_ROTATION_AZIMUTH']

mesh_central_lat = float(mesh_central_lat.lower().replace('d', 'e'))
mesh_central_lon = float(mesh_central_lon.lower().replace('d', 'e'))
mesh_width_xi = float(mesh_width_xi.lower().replace('d', 'e'))
mesh_width_eta = float(mesh_width_eta.lower().replace('d', 'e'))
mesh_gamma_rot = float(mesh_gamma_rot.lower().replace('d', 'e'))



df = pd.read_csv('events_with_xi_eta.csv')
df['time'] = pd.to_datetime(df['time'])

width_xi = 51
width_eta = 51
dxi = 1
deta = 1

n_xi = width_xi // dxi + 1
n_eta = width_eta // deta + 1
dxi = width_xi / n_xi
deta = width_eta / n_eta

grids_xi0 = -width_xi/2 + np.arange(n_xi) * dxi
grids_eta0 = -width_eta/2 + np.arange(n_eta) * deta

for xi0 in grids_xi0:
    for eta0 in grids_eta0:

        xi1 = xi0 + dxi
        eta1 = eta0 + deta

        mask = (df['xi'] >= xi0) & (df['xi'] < xi1) & (df['eta'] >= eta0) & (df['eta'] < eta1)
        if not mask.any():
            continue
        df_sel = df[mask]

        df_sel = df_sel.sort_values(by=['mag', 'time'], ascending=False)

        # print(f"{xi0=}, {xi1=}, {eta0=}, {eta1=}")
        # print("==========================")
        # print(df_sel)
        print(df_sel.iloc[0]['name'])

        # mask = df_sel['mag'] < 6.5
        # if any(mask):
        #     print(df_sel[mask].iloc[0]['name'])
