import os
import pandas as pd
import numpy as np
import pyvista as pv
import argparse

from meshfem3d_utils import xyz2latlon_deg, geodetic_lat2geocentric_lat

parser = argparse.ArgumentParser()

parser.add_argument("par_file", help="SEM Par_file")
parser.add_argument(
    "-n",
    "--ngrid",
    nargs=3,
    metavar=("xi", "eta", "r"),
    help="number of grids along %(metavar)s",
    default=[11, 11, 11],
    type=int,
)
parser.add_argument(
    "-r",
    "--r_range",
    nargs=2,
    metavar=("min_r", "max_r"),
    default=[0.8, 1.0],
    type=float,
)
parser.add_argument("--vtk", default="pointspread_points.vtk", help="output VTK file")

parser.add_argument("-o", "--out_dir", default="output", help="output directory")

parser.add_argument(
    "--vertical_list", default="vertical_slices.csv", help="filename of vertical slices list"
)

parser.add_argument(
    "--horizontal_list", default="horizontal_slices.csv", help="filename of vertical slices list"
)

args = parser.parse_args()
print(args)

sem_parfile = args.par_file
vtk_file = f"{args.out_dir}/{args.vtk}"
# slice_file = args.slice

sem_params = pd.read_csv(
    sem_parfile,
    delimiter=r"\s*=\s*",
    header=None,
    comment="#",
    names=["key", "value"],
    dtype=dict(key=object, value=object),
    index_col=["key"],
    engine="python",
).to_dict()["value"]

mesh_central_lat = sem_params["CENTER_LATITUDE_IN_DEGREES"]
mesh_central_lon = sem_params["CENTER_LONGITUDE_IN_DEGREES"]
mesh_width_xi = sem_params["ANGULAR_WIDTH_XI_IN_DEGREES"]
mesh_width_eta = sem_params["ANGULAR_WIDTH_ETA_IN_DEGREES"]
mesh_gamma_rot = sem_params["GAMMA_ROTATION_AZIMUTH"]

lat_center = float(mesh_central_lat.lower().replace("d", "e"))
lon_center = float(mesh_central_lon.lower().replace("d", "e"))
xi_width = float(mesh_width_xi.lower().replace("d", "e"))
eta_width = float(mesh_width_eta.lower().replace("d", "e"))
gamma_rot = float(mesh_gamma_rot.lower().replace("d", "e"))

lat0 = np.deg2rad(lat_center)
lon0 = np.deg2rad(lon_center)

theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat0)
# print(lat0, theta)
phi = lon0

v0_r = np.array(
    [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
)
v0_e = np.array([-np.sin(phi), np.cos(phi), 0])
v0_n = np.array(
    [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
)

# rotate (v_east, v_north) thourgh v_radial by gamma_rot to (v_xi, v_eta)
gamma = np.deg2rad(gamma_rot)
v0_xi = np.cos(gamma) * v0_e + np.sin(gamma) * v0_n
v0_eta = -np.sin(gamma) * v0_e + np.cos(gamma) * v0_n

# points
n_xi, n_eta, n_r = args.ngrid
npts = n_xi * n_eta * n_r
points = np.zeros((n_xi, n_eta, n_r, 3))

rmin, rmax = args.r_range
r_grid = np.linspace(rmin, rmax, n_r)

angle_xi = np.deg2rad(xi_width)
angle_eta = np.deg2rad(eta_width)

half_angle_xi = 0.5 * angle_xi
half_angle_eta = 0.5 * angle_eta

dxi = angle_xi / (n_xi - 1)
deta = angle_eta / (n_eta - 1)

for ixi in np.arange(n_xi):
    xi = -half_angle_xi + dxi * ixi
    for ieta in np.arange(n_eta):
        eta = -half_angle_eta + deta * ieta

        l_xi = np.tan(xi)
        l_eta = np.tan(eta)

        v = v0_xi * l_xi + v0_eta * l_eta + v0_r
        v = v / sum(v**2) ** 0.5

        points[ixi, ieta, :, :] = r_grid[:, None] * v[None, :]


mesh = pv.PolyData(points.reshape((-1, 3)))
mesh.save(vtk_file)

# make cross-sections
slice_params = []

# vertical xsections (great circle plane) parallel to v0_xi
for ieta in np.arange(n_eta):
    xi = 0
    eta = -half_angle_eta + deta * ieta

    r = v0_eta * np.sin(eta) + v0_r * np.cos(eta)
    r = r / sum(r**2) ** 0.5
    lat, lon = xyz2latlon_deg(r)

    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat)
    phi = lon
    e = np.array([-np.sin(phi), np.cos(phi), 0])
    n = np.array(
        [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
    )

    half_angle_theta = np.arctan(np.cos(eta) * np.tan(half_angle_xi))

    azimuth = np.arctan2(np.dot(v0_xi, e), np.dot(v0_xi, n))

    params = (
        np.rad2deg(xi),
        np.rad2deg(eta),
        np.rad2deg(lat),
        np.rad2deg(lon),
        np.rad2deg(azimuth),
        -np.rad2deg(half_angle_theta),
        np.rad2deg(half_angle_theta),
    )
    slice_params.append(params)

# vertical xsections (great circle plane) parallel to v0_eta
for ixi in np.arange(n_xi):
    xi = -half_angle_xi + dxi * ixi
    eta = 0

    r = v0_xi * np.sin(xi) + v0_r * np.cos(xi)
    r = r / sum(r**2) ** 0.5
    lat, lon = xyz2latlon_deg(r)

    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat)
    phi = lon
    e = np.array([-np.sin(phi), np.cos(phi), 0])
    n = np.array(
        [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
    )

    half_angle_theta = np.arctan(np.cos(xi) * np.tan(half_angle_eta))

    azimuth = np.arctan2(np.dot(v0_eta, e), np.dot(v0_eta, n))

    params = (
        np.rad2deg(xi),
        np.rad2deg(eta),
        np.rad2deg(lat),
        np.rad2deg(lon),
        np.rad2deg(azimuth),
        -np.rad2deg(half_angle_theta),
        np.rad2deg(half_angle_theta),
    )
    slice_params.append(params)

column_names = ["xi", "eta", "lat", "lon", "azimuth", "min_theta", "max_theta"]
df = pd.DataFrame(slice_params, columns=column_names)
df.to_csv(os.path.join(args.out_dir, args.vertical_list), float_format="%15.5e", index=False)


# make horizontal cross-sections (spherical cap)
column_names = ["radius", "central_lat", "central_lon", "width_xi", "width_eta", "rotation_angle"]
slice_params = [ (r, lat_center, lon_center, xi_width, eta_width, gamma_rot) for r in r_grid]
df = pd.DataFrame(slice_params, columns=column_names)
df.to_csv(os.path.join(args.out_dir, args.horizontal_list), float_format="%15.5e", index=False)

# with open(slice_file, "w") as fp:
#     fp.write("xi  eta  lat  lon  azimuth  min_theta max_theta\n")
#     for params in slice_params:
#         out = "%8.3f %8.3f %10.3e %10.3e %10.3e %8.3f %8.3f\n" % (params)
#         fp.write(out)
