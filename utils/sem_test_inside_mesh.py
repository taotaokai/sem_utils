import sys
import numpy as np
import pandas as pd
from misfit.misfit import get_dist_from_mesh_boundary

def vector2latlon_deg(v):
    lon = np.arctan2(v[1], v[0])
    latc = np.arctan2(v[2], (v[0]**2+v[1]**2)**0.5)
    f = 1/299.8
    factor = 1.0 / (1 - f)**2
    lat = np.arctan(factor * np.tan(latc))
    return np.rad2deg(lat), np.rad2deg(lon)

def geodetic_lat2geocentric_lat(geodetic_lat):
    f = 1/299.8
    factor = (1 - f)**2
    return np.arctan(factor * np.tan(geodetic_lat))

def check_bound(lat_center, lon_center, gamma_rot, lat_test, lon_test):
    lat0 = np.deg2rad(lat_center)
    lon0 = np.deg2rad(lon_center)
    # assume zero altitude from the reference ellipsoid
    theta = 0.5*np.pi - geodetic_lat2geocentric_lat(lat0)
    phi = lon0
    # radial/easting/northing direction at (lat_center, lon_center)
    v0_r = np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])
    v0_e = np.array([-np.sin(phi), np.cos(phi), 0])
    v0_n = np.array([-np.cos(theta) * np.cos(phi),
                     -np.cos(theta) * np.sin(phi),
                     np.sin(theta)])

    # rotate (v0_e, v0_n) to (v0_xi, v0_eta) through v0_r by gamma_rot counter-clockwise
    gamma = np.deg2rad(gamma_rot)
    v0_xi = np.cos(gamma) * v0_e + np.sin(gamma) * v0_n
    v0_eta = -np.sin(gamma) * v0_e + np.cos(gamma) * v0_n
    # print(v_xi, v_eta)

    # test point
    lat1 = np.deg2rad(lat_test)
    lon1 = np.deg2rad(lon_test)
    theta = 0.5*np.pi - geodetic_lat2geocentric_lat(lat1)
    phi = lon1
    v1 = np.array([np.sin(theta) * np.cos(phi),
                   np.sin(theta) * np.sin(phi),
                   np.cos(theta)])
    # project to v0_r, v0_xi, v0_eta
    l_r = np.dot(v0_r, v1)
    l_xi = np.dot(v0_xi, v1)
    l_eta = np.dot(v0_eta, v1)
    angle_xi = np.rad2deg(np.arctan2(l_xi, l_r))
    angle_eta = np.rad2deg(np.arctan2(l_eta, l_r))

    return angle_xi, angle_eta


if __name__ == '__main__':
    # Usage: cat LIST | ./sem_test_inside_mesh Par_file
    # LIST is a file containing "lat lon name" on each row (like, 135 35 ST1).
    # each row of output: "name min_distance_from_mesh_boundary (negative: outside)"

    import argparse
    # 1. Create a parser
    parser = argparse.ArgumentParser(description='A simple script demonstrating argparse.')
    # 2. Add arguments
    parser.add_argument('sem_parfile', type=str, help='SEM Par_file.')
    # parser.add_argument('--output', '-o', type=str, default='output.txt',
    #                     help='The path to the output file (default: output.txt).')
    # parser.add_argument('--verbose', '-v', action='store_true',
    #                     help='Enable verbose output.')
    # 3. Parse arguments
    args = parser.parse_args()

    sem_parfile = args.sem_parfile

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
        xi, eta = get_dist_from_mesh_boundary(mesh_central_lat, mesh_central_lon, mesh_gamma_rot, float(lat), float(lon))
        dist_xi = 0.5 * mesh_width_xi - abs(xi)
        dist_eta = 0.5 * mesh_width_eta - abs(eta)
        min_dist = min(dist_xi, dist_eta)
        print(f"{name} {min_dist:6.2f}")
