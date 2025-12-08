import sys
import pandas as pd
import numpy as np
import pyvista as pv
# import simplekml
import argparse

from meshfem3d_utils import geodetic_lat2geocentric_lat

parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int, help="number of mesh mpi processes")
parser.add_argument("xsection_list", help="file of xsection params")
parser.add_argument(
    "-d",
    "--sem_db",
    help="SEM DATABASES_MPI",
    default="DATABASES_MPI",
)
parser.add_argument(
    "-m",
    "--model",
    nargs="+",
    default=["vpv"],
)
parser.add_argument("-o", "--out", default="slices.nc", help="output NETCDF file")
parser.add_argument(
    "-n",
    "--ngrid",
    nargs=2,
    metavar=("xi", "eta"),
    help="number of grids along %(metavar)s",
    default=[100, 100],
    type=int,
)
parser.add_argument("--vtk", action='store_true', help='Output VTK file')

args = parser.parse_args()
print(args)

xsection_params = pd.read_csv(args.xsection_list)

nxi, neta = args.ngrid

for i, params in xsection_params.iterrows():

    radius = params['radius']
    lat = np.deg2rad(params['central_lat'])
    lon = np.deg2rad(params['central_lon'])
    gamma = np.deg2rad(params["rotation_angle"])
    angle_xi = np.deg2rad(params["width_xi"]) 
    angle_eta = np.deg2rad(params["width_eta"]) 

    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat)
    phi = lon
    ve = np.array([-np.sin(phi), np.cos(phi), 0])
    vn = np.array(
        [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
    )
    vr = np.array(
        [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )

    v_xi = np.cos(gamma) * ve + np.sin(gamma) * vn
    v_eta = -np.sin(gamma) * ve + np.cos(gamma) * vn

    half_angle_xi = 0.5 * angle_xi
    half_angle_eta = 0.5 * angle_eta

    dxi = angle_xi / (nxi - 1)
    deta = angle_eta / (neta - 1)

    xi = -half_angle_xi + dxi * np.arange(nxi)
    eta = -half_angle_eta + deta * np.arange(neta)
    xi2, eta2 = np.meshgrid(xi, eta, indexing='ij')

    l_xi = np.tan(xi2)
    l_eta = np.tan(eta2)
    v = l_xi[:,:,None] * v_xi[None,None,:] + l_eta[:,:,None] * v_eta[None,None,:] + vr[None,None,:]
    v = v / np.sum(v**2,axis=-1,keepdims=True)** 0.5
    points = radius * v

    if args.vtk:
        # connectivity
        ncells = (nxi - 1) * (neta - 1)
        connectivity = np.zeros((ncells, 4), dtype=int)
        iy, ix = np.unravel_index(np.arange(ncells), (nxi - 1, neta - 1))
        ii = 0
        for dx, dy in (0, 0), (1, 0), (1, 1), (0, 1):
            ind = np.ravel_multi_index((iy + dy, ix + dx), (nxi, neta))
            connectivity[:, ii] = (ix + dx) + neta * (iy + dy)
            ii += 1
        mesh = pv.UnstructuredGrid(
            {pv.CellType.QUAD: connectivity}, points.reshape((-1, 3))
        )
        mesh.save(f"horz_xsection_{i:02d}.vtk")