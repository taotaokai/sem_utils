import sys
import numpy as np
import pandas as pd

from meshfem3d_utils import sem_latlon2xieta

# def vector2latlon_deg(v):
#     lon = np.arctan2(v[1], v[0])
#     latc = np.arctan2(v[2], (v[0]**2+v[1]**2)**0.5)
#     f = 1/299.8
#     factor = 1.0 / (1 - f)**2
#     lat = np.arctan(factor * np.tan(latc))
#     return np.rad2deg(lat), np.rad2deg(lon)

if __name__ == '__main__':

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

    for line in sys.stdin:
        line = line.strip()
        lat, lon, name = line.strip().split()
        xi, eta = sem_latlon2xieta(mesh_central_lat, mesh_central_lon, mesh_gamma_rot, float(lat), float(lon))
        dist_xi = 0.5 * mesh_width_xi - abs(xi)
        dist_eta = 0.5 * mesh_width_eta - abs(eta)
        min_dist = min(dist_xi, dist_eta)
        print(f"{name:15s} {min_dist:6.2f} {xi:6.2f} {eta:6.2f}")
