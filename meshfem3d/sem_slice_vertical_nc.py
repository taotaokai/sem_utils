import sys
import pandas as pd
import numpy as np
import pyvista as pv
import pyproj

# import simplekml
import argparse
import xarray as xr
from netCDF4 import Dataset
from scipy.interpolate import interp1d, interpn, RegularGridInterpolator

from meshfem3d_utils import xyz2latlon_deg, geodetic_lat2geocentric_lat
from meshfem3d_constants import R_EARTH

parser = argparse.ArgumentParser()

# parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
parser.add_argument("xsection_list", help="file of xsection params")
parser.add_argument(
    "-f",
    "--model_file",
    help="model nc file",
)
parser.add_argument(
    "-m",
    "--model_name",
    # nargs="+",
    # default=["vpv"],
    default="VPV",
)
parser.add_argument("-o", "--out_dir", default="output", help="output directory")
parser.add_argument(
    "-r",
    "--r_range",
    nargs=2,
    metavar=("rmin", "rmax"),
    help="limits along radial direction",
    default=[0.8, 1.0],
    type=float,
)
parser.add_argument(
    "-n",
    "--ngrid",
    nargs=2,
    metavar=("theta", "r"),
    help="number of grids along %(metavar)s",
    default=[100, 50],
    type=int,
)
parser.add_argument("--vtk", action="store_true", help="Output VTK file")

args = parser.parse_args()
print(args)

ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")  # ECEF to GPS

#====== read in topo file
bedrock = Dataset("ETOPO/ETOPO_2022_v1_60s_N90W180_bed.nc", "r", format="NETCDF4")
geoid = Dataset("ETOPO/ETOPO_2022_v1_60s_N90W180_geoid.nc", "r", format="NETCDF4")
etopo1_heights = np.array(bedrock.variables['z']) + np.array(geoid.variables['z']) # z(lat, lon) height from WGS84 ellipsoid
etopo1_lons = np.array(geoid.variables['lon'])
etopo1_lats = np.array(geoid.variables['lat'])
etopo1_interp = RegularGridInterpolator((etopo1_lats, etopo1_lons), etopo1_heights, bounds_error=False, fill_value=np.nan)
bedrock.close()
geoid.close()

#====== read in model file
model = xr.open_dataset(args.model_file)
model_depths = np.array(model.variables['depth'])
model_lats = np.array(model.variables['latitude'])
model_lons = np.array(model.variables['longitude'])
model_values = np.array(model.variables[args.model_name]) / 1000.0 # m/s to km/s
nlat, nlon = model_lats.size, model_lons.size
model_values = model_values.reshape((-1, nlon, nlat))
model_interp = RegularGridInterpolator((model_depths, model_lons, model_lats), model_values, bounds_error=False, fill_value=np.nan)
model.close()

#====== read in xsection file
xsection_params = pd.read_csv(args.xsection_list)

ntheta, nr = args.ngrid
rmin, rmax = args.r_range
r_grid = np.linspace(rmin, rmax, nr)



for i, params in xsection_params.iterrows():

    lat = np.deg2rad(params["lat"])
    lon = np.deg2rad(params["lon"])
    azimuth = np.deg2rad(params["azimuth"])
    min_theta = np.deg2rad(params["min_theta"])
    max_theta = np.deg2rad(params["max_theta"])

    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat)
    phi = lon
    ve = np.array([-np.sin(phi), np.cos(phi), 0])
    vn = np.array(
        [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
    )
    vr = np.array(
        [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    vx = np.cos(azimuth) * vn + np.sin(azimuth) * ve  # vector along xsection

    theta_width = max_theta - min_theta
    dtheta = theta_width / (ntheta - 1)
    theta_grid = min_theta + dtheta * np.arange(ntheta)

    points = np.zeros((ntheta, nr, 3))
    tangents = np.zeros((ntheta, 3))
    for itheta, theta in enumerate(theta_grid):
        # theta = min_theta + itheta * dtheta
        v = np.tan(theta) * vx + vr
        v = v / sum(v**2) ** 0.5
        t = vx - sum(v * vx) * v
        t = t / sum(t**2) ** 0.5
        points[itheta, :, :] = r_grid[:, None] * v[None, :]
        tangents[itheta, :] = t

    #----- interpolation
    # convert to lat,lon,alt
    xyz = R_EARTH * points.reshape((-1, 3))
    lat, lon, alt = ecef2gps.transform(xyz[:,0], xyz[:,1], xyz[:,2])
    print(np.min(alt), np.max(alt))
    print(np.min(lat), np.max(lat))
    print(np.min(lon), np.max(lon))

    topo = etopo1_interp((lat, lon))
    print(np.min(topo), np.max(topo))

    # convert alt(ellipsoidal height) to depth_km below surface
    depth = -1 * (alt - topo) / 1000.0 # to km

    print(np.min(depth), np.max(depth))

    # depth[depth < 0.001] = 0.001

    # interpolate model values
    output_model = model_interp((depth, lon, lat))
    output_model = output_model.reshape((ntheta, nr))

    if args.vtk:
        blocks = pv.MultiBlock()
        # connectivity
        ncells = (ntheta - 1) * (nr - 1)
        connectivity = np.zeros((ncells, 4), dtype=int)
        iy, ix = np.unravel_index(np.arange(ncells), (ntheta - 1, nr - 1))
        ii = 0
        for dx, dy in (0, 0), (1, 0), (1, 1), (0, 1):
            ind = np.ravel_multi_index((iy + dy, ix + dx), (ntheta, nr))
            connectivity[:, ii] = (ix + dx) + nr * (iy + dy)
            ii += 1
        mesh = pv.UnstructuredGrid(
            {pv.CellType.QUAD: connectivity}, points.reshape((-1, 3))
        )
        mesh.point_data[args.model_name] = output_model.flatten()
        blocks.append(mesh)
        mean_val = np.nanmean(output_model)
        scale = 2 * 10 * dtheta / mean_val
        # plot radial 1-D profile
        for itheta in range(10, ntheta, 10):
            m0 = output_model[itheta, :]
            mask = np.isnan(m0)
            if np.all(mask): # skip profile of all nan's
                continue
            m = output_model[itheta, :] - mean_val
            p = points[itheta, :, :]
            t = tangents[itheta, :]
            x = p + scale * m[:, None] * t[None, :]
            line = pv.lines_from_points(x)
            # line.point_data[args.model_name] = m0
            blocks.append(line)
            line = pv.lines_from_points(p)
            line.point_data[args.model_name] = m0
            blocks.append(line)
        blocks.save(f"{args.out_dir}/vert_xsection_{i:02d}.vtmb")

    #
    # ds_xsection = xr.Dataset()
    # da = xr.DataArray(
    #     output_model,
    #     dims=("theta", "radius"),
    #     coords={
    #         "theta": ("theta", np.rad2deg(theta_grid), {"unit": "degree"}),
    #         "radius": (
    #             "radius",
    #             r_grid,
    #             {"unit": "R_EARTH"},
    #         ),
    #         # "dim3": ("dim3", [0,1,2], {"unit":"xyz"}),
    #     },
    #     attrs={"latc": params["lat"], "lonc": params["lon"], "azimuth": params["azimuth"]},
    #     name=args.model_name,
    # )
    # da.to_netcdf(f"{args.out_dir}/vert_xsection_{i:02d}.nc")
    # ds_xsection[f"vertical_xsection_{i:02d}"] = da

# ds_xsection.to_netcdf("vertical_xsections.netcdf")
