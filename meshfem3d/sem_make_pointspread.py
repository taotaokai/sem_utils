import os
import pandas as pd
import numpy as np
import pyvista as pv
import argparse

from meshfem3d_utils import R_EARTH_KM, ecef2latlon_zeroalt, geodetic_lat2geocentric_lat

parser = argparse.ArgumentParser()

parser.add_argument("center_lat", help="mesh center's latitude in degrees", type=float)
parser.add_argument("center_lon", help="mesh center's longitude in degrees", type=float)
parser.add_argument(
    "width_xi",
    help="width of xi in degrees, easting if rotation angle is zero",
    type=float,
)
parser.add_argument(
    "width_eta",
    help="width of eta in degrees, northing if rotation angle is zero",
    type=float,
)
parser.add_argument(
    "rotation_angle", help="rotation angle in degrees (anti-clockwise)", type=float
)

# parser.add_argument("--xi_range", nargs=2, default=[-1, 1], help="range of xi in degrees, easting if rotation angle is zero", type=float)
# parser.add_argument("--n_xi", help="number of grids along xi", type=int)
# parser.add_argument("--eta_range", nargs=2, default=[-1, 1], help="range of eta in degrees, easting if rotation angle is zero", type=float)
# parser.add_argument("--n_eta", help="number of grids along eta", type=int)
# parser.add_argument("--depth_range", nargs=2, default=[10, 20], help="range of depth in km", type=float)
# parser.add_argument("--n_depth", help="number of grids along xi", type=int)

parser.add_argument(
    "--xi",
    nargs="+",
    default=[-1, 1],
    help="grids of xi in degrees, easting if rotation angle is zero",
    type=float,
)
parser.add_argument(
    "--eta",
    nargs="+",
    default=[-1, 1],
    help="grids of eta in degrees, northing if rotation angle is zero",
    type=float,
)
parser.add_argument(
    "--depth", nargs="+", default=[30, 100], help="grids of depth in km", type=float
)

parser.add_argument("--out_dir", default="./", help="output directory")
parser.add_argument("--vtk", default="pointspread_points.vtk", help="output VTK file")
parser.add_argument(
    "--point_list",
    default="pointspread_points.csv",
    help="filename of point-spread point list",
)
parser.add_argument(
    "--vertical_list",
    default="vertical_slices.csv",
    help="filename of vertical slices list",
)
parser.add_argument(
    "--horizontal_list",
    default="horizontal_slices.csv",
    help="filename of vertical slices list",
)

args = parser.parse_args()
print(args)

# mesh center
lat0 = np.deg2rad(args.center_lat)
lon0 = np.deg2rad(args.center_lon)

theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat0)
phi = lon0

v0_r = np.array(
    [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
)
v0_e = np.array([-np.sin(phi), np.cos(phi), 0])
v0_n = np.array(
    [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
)

# rotate (v_east, v_north) thourgh v_radial by gamma_rot to (v_xi, v_eta)
gamma = np.deg2rad(args.rotation_angle)
v0_xi = np.cos(gamma) * v0_e + np.sin(gamma) * v0_n
v0_eta = -np.sin(gamma) * v0_e + np.cos(gamma) * v0_n

# grid locations for point-spread function
grid_xi = np.deg2rad(args.xi)
grid_eta = np.deg2rad(args.eta)
grid_depth = np.array(args.depth)
grid_radius = 1.0 - grid_depth / R_EARTH_KM

n_xi, n_eta, n_depth = len(grid_xi), len(grid_eta), len(grid_depth)
npts = n_xi * n_eta * n_depth
points = np.zeros((n_xi, n_eta, n_depth, 3))

# grid_xi = np.deg2rad(np.linspace(args.xi_range[0], args.xi_range[1], n_xi))
# grid_eta = np.deg2rad(np.linspace(args.eta_range[0], args.eta_range[1], n_eta))
# grid_depth = np.linspace(args.depth_range[0], args.depth_range[1], n_depth)
# grid_radius = 1.0 - grid_depth / R_EARTH_KM

data = []
for i, xi in enumerate(grid_xi):
    for j, eta in enumerate(grid_eta):
        l_xi = np.tan(xi)
        l_eta = np.tan(eta)
        v = v0_xi * l_xi + v0_eta * l_eta + v0_r
        v = v / sum(v**2) ** 0.5
        for k, depth in enumerate(grid_depth):
            points[i, j, k, :] = grid_radius[k] * v
            data.append(
                {
                    "xi": args.xi[i],
                    "eta": args.eta[j],
                    "depth": args.depth[k],
                    "x": points[i, j, k, 0],
                    "y": points[i, j, k, 1],
                    "z": points[i, j, k, 2],
                }
            )
        # points[i, j, :, :] = grid_radius[:, None] * v[None, :]

# write to csv
df = pd.DataFrame(data)
df.to_csv(
    os.path.join(args.out_dir, args.point_list), float_format="%15.5e", index=False
)

# save to vtk
mesh = pv.PolyData(points.reshape((-1, 3)))
mesh.save(os.path.join(args.out_dir, args.vtk))

# make cross-sections
slice_params = []

half_angle_xi = np.deg2rad(args.width_xi / 2)
half_angle_eta = np.deg2rad(args.width_eta / 2)

# vertical xsections (great circle plane) parallel to v0_xi
for eta in grid_eta:
    xi = 0

    r = v0_eta * np.sin(eta) + v0_r * np.cos(eta)
    r = r / sum(r**2) ** 0.5
    lat, lon = ecef2latlon_zeroalt(r[0], r[1], r[2])

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
for xi in grid_xi:
    eta = 0

    r = v0_xi * np.sin(xi) + v0_r * np.cos(xi)
    r = r / sum(r**2) ** 0.5
    lat, lon = ecef2latlon_zeroalt(r[0], r[1], r[2])

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
df.to_csv(
    os.path.join(args.out_dir, args.vertical_list), float_format="%15.5e", index=False
)


# make horizontal cross-sections (spherical cap)
column_names = [
    "depth_km",
    "central_lat",
    "central_lon",
    "width_xi",
    "width_eta",
    "rotation_angle",
]
slice_params = [
    (
        depth,
        args.center_lat,
        args.center_lon,
        args.width_xi,
        args.width_eta,
        args.rotation_angle,
    )
    for depth in grid_depth
]
df = pd.DataFrame(slice_params, columns=column_names)
df.to_csv(
    os.path.join(args.out_dir, args.horizontal_list), float_format="%15.5e", index=False
)

# with open(slice_file, "w") as fp:
#     fp.write("xi  eta  lat  lon  azimuth  min_theta max_theta\n")
#     for params in slice_params:
#         out = "%8.3f %8.3f %10.3e %10.3e %10.3e %8.3f %8.3f\n" % (params)
#         fp.write(out)
