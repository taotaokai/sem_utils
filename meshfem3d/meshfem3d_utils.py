#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import warnings

import numpy as np

from meshfem3d_constants import *

import numba

# ///////////////////////////////////////////////
# constants.h
# NGLLX = 5
# NGLLY = NGLLX
# NGLLZ = NGLLX
#
# MIDX = int((NGLLX-1)/2)
# MIDY = int((NGLLY-1)/2)
# MIDZ = int((NGLLZ-1)/2)
#
# GAUSSALPHA = 0
# GAUSSBETA = 0
#
## proc000*_reg1_solver_data.bin
# MESH_ARRAY_LIST = [
#  ('nspec','i4')                      ,
#  ('nglob','i4')                      ,
#  ('x','f4')                          ,
#  ('y','f4')                          ,
#  ('z','f4')                          ,
#  ('ibool','i4')                      ,
#  ('idoubling','i4')                  ,
#  ('ispec_is_tiso','i4')              ,
#  ('dxsi_dx','f4')                      ,
#  ('dxsi_dy','f4')                      ,
#  ('dxsi_dz','f4')                      ,
#  ('deta_dx','f4')                     ,
#  ('deta_dy','f4')                     ,
#  ('deta_dz','f4')                     ,
#  ('dgam_dx','f4')                   ,
#  ('dgam_dy','f4')                   ,
#  ('dgam_dz','f4')                   ,
#  ]

#==================================================#

@numba.jit("Tuple((float64[:], float64[:]))(int64)", nogil=True) #, cache=True)
def gll_nodes_weights(n :int):
    if n == 5:
        x = np.array([-1, -((3.0 / 7.0) ** 0.5), 0, (3.0 / 7.0) ** 0.5, 1])
        w = np.array([0.1, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 0.1])
        return x, w
    else:
        # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
        x = -1 * np.cos(np.pi * np.arange(n) / (n-1))
        # n1 = n - 1
        # The Legendre Vandermonde Matrix
        P = np.zeros((n, n))
        # Compute P_(N) using the recursion relation
        # Compute its first and second derivatives and
        # update x using the Newton-Raphson method.
        xold = 2 * np.ones_like(x)
        while np.max(np.abs(x - xold)) > 1e-7:
            xold = x
            P[:, 0] = 1
            P[:, 1] = x
            for k in range(1, n-1):
                P[:, k + 1] = ((2 * k + 1) * x * P[:, k] - (k) * P[:, k - 1]) / (k + 1)
            x = xold - (x * P[:, -1] - P[:, -2]) / (n * P[:, -1])
        w = 2 / (n * (n -1) * P[:, -1] ** 2)
        return x, w

@numba.jit("float64[:](float64[:], float64)", nogil=True) #, cache=True)
def lagrange_poly(nodes, x:float):
    """
    n-node Lagrange basis evaluated at x
    """
    z = nodes
    n = len(nodes)
    # lag_i(x)
    lag = np.zeros(n)
    ind = np.arange(n)
    for i in range(n):
        zj = z[ind != i]
        lag[i] = np.prod(x - zj) / np.prod(z[i] - zj)
    return lag

@numba.jit("float64[:,:](float64[:], float64[:])", nogil=True) #, cache=True)
def lagrange_poly_derivative(nodes, xx):
    """
    derivative of n-node Lagrange basis at xx
    """
    z = nodes
    nz = len(nodes)
    nx = len(xx)
    # a[i,j] = dlag_i/dx(x_j)
    dlag_dx = np.zeros((nz, nx))
    ind = np.arange(nz)
    for i in range(nz):
        denominator = np.prod((z[i] - z[ind != i]))
        for j in range(nx):
            x = xx[j]
            numerator = 0
            for k in range(nz):
                if k == i:
                    continue
                mask = (ind != i) & (ind != k)
                numerator += np.prod(x - z[mask])
            dlag_dx[i, j] = numerator / denominator
    return dlag_dx

def get_gll_weights():
    """ use sem build-in gll_library
    """
    import gll_library
    # GLL points and weights
    zgll, wgll = gll_library.zwgljd(NGLLX, GAUSSALPHA, GAUSSBETA)
    # Get derivatives of langrange basis at gll nodes
    dlag_dzgll = np.zeros((NGLLX, NGLLX))
    for i in range(NGLLX):
        for j in range(NGLLX):
            # dlagP_i/dx(xgll[j])
            dlag_dzgll[i, j] = gll_library.lagrange_deriv_gll(i, j, zgll, NGLLX)
    return zgll, wgll, dlag_dzgll

#==================================================#

def geodetic_lat2geocentric_lat(geodetic_lat):
    f = EARTH_FLATTENING
    factor = (1 - f) ** 2
    return np.arctan(factor * np.tan(geodetic_lat))


def xyz2latlon_deg(v):
    lon = np.arctan2(v[1], v[0])
    latc = np.arctan2(v[2], (v[0] ** 2 + v[1] ** 2) ** 0.5)
    f = EARTH_FLATTENING
    factor = 1.0 / (1 - f) ** 2
    lat = np.arctan(factor * np.tan(latc))
    return lat, lon


def sem_latlon2xieta(lat_center, lon_center, gamma_rot, lat_test, lon_test):
    """convert lat/lon to SEM mesh local coordinates xi/eta"""
    lat0 = np.deg2rad(lat_center)
    lon0 = np.deg2rad(lon_center)
    # assume zero altitude from the reference ellipsoid
    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat0)
    phi = lon0
    # radial/easting/northing direction at (lat_center, lon_center)
    v0_r = np.array(
        [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    v0_e = np.array([-np.sin(phi), np.cos(phi), 0])
    v0_n = np.array(
        [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
    )

    # rotate (v0_e, v0_n) to (v0_xi, v0_eta) through v0_r by gamma_rot counter-clockwise
    gamma = np.deg2rad(gamma_rot)
    v0_xi = np.cos(gamma) * v0_e + np.sin(gamma) * v0_n
    v0_eta = -np.sin(gamma) * v0_e + np.cos(gamma) * v0_n
    # print(v_xi, v_eta)

    # conver test point to (xi, eta)
    lat1 = np.deg2rad(lat_test)
    lon1 = np.deg2rad(lon_test)
    theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat1)
    phi = lon1
    v1 = np.array(
        [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    # project to v0_r, v0_xi, v0_eta
    l_r = np.dot(v0_r, v1)
    l_xi = np.dot(v0_xi, v1)
    l_eta = np.dot(v0_eta, v1)
    angle_xi = np.rad2deg(np.arctan2(l_xi, l_r))
    angle_eta = np.rad2deg(np.arctan2(l_eta, l_r))

    return angle_xi, angle_eta


def rotmat_enu_to_ecef(lon, lat):
    """rotation matrix from local ENU to ECEF coordinate basises
    rotmat[:,0] = Ve # column vector of the Easting direction in ECEF coordinate
    rotmat[:,1] = Vn # column vector of the Northing direction in ECEF coordinate
    rotmat[:,2] = Vu # column vector of the Up (ellipsoid height) direction in ECEF coordinate

    xyz_ecef = xyz0_ecef + rotmat * enu
    enu = transpose(rotmat) * (xyz_ecef - xyz0_ecef)

    , where xyz0_ecef is the reference point at (lon,lat,alt).
    """
    coslat = np.cos(np.deg2rad(lat))
    sinlat = np.sin(np.deg2rad(lat))
    coslon = np.cos(np.deg2rad(lon))
    sinlon = np.sin(np.deg2rad(lon))

    rotmat = np.zeros((3, 3))
    rotmat[0, :] = [-sinlon, -sinlat * coslon, coslat * coslon]
    rotmat[1, :] = [coslon, -sinlat * sinlon, coslat * sinlon]
    rotmat[2, :] = [0.0, coslat, sinlat]

    return rotmat


def sem_mesh_read(mesh_file):
    """read in SEM solver_data.bin slice"""
    from scipy.io import FortranFile

    mesh_data = {}

    with FortranFile(mesh_file, "r") as f:
        for field in MESH_ARRAY_LIST:
            field_name = field[0]
            data_type = field[1]
            mesh_data[field_name] = f.read_ints(dtype=data_type)

    mesh_data["nspec"] = mesh_data["nspec"][0]
    mesh_data["nglob"] = mesh_data["nglob"][0]

    # GLL dims
    gll_dims = (mesh_data["nspec"], NGLLZ, NGLLY, NGLLX)
    mesh_data["gll_dims"] = gll_dims

    # reshape
    for field_name in [
        "ibool",
        "dxsi_dx",
        "dxsi_dy",
        "dxsi_dz",
        "deta_dx",
        "deta_dy",
        "deta_dz",
        "dgam_dx",
        "dgam_dy",
        "dgam_dz",
    ]:
        # NB: binary files are written in Fortran column-major convention !!!
        # NB: reshape 1-D array to 4-D GLL tensor by Fortran convention, and
        # NB: also convert to a contiguous array in memory, in case of direct memory copy or MPI transfer
        mesh_data[field_name] = np.reshape(mesh_data[field_name], gll_dims)

    # jacobian: det( d(x,y,z)/d(xi,eta,gamma))
    mesh_data["jacobian"] = 1.0 / (
        mesh_data["dxsi_dx"]
        * (
            mesh_data["deta_dy"] * mesh_data["dgam_dz"]
            - mesh_data["deta_dz"] * mesh_data["dgam_dy"]
        )
        - mesh_data["dxsi_dy"]
        * (
            mesh_data["deta_dx"] * mesh_data["dgam_dz"]
            - mesh_data["deta_dz"] * mesh_data["dgam_dx"]
        )
        + mesh_data["dxsi_dz"]
        * (
            mesh_data["deta_dx"] * mesh_data["dgam_dy"]
            - mesh_data["deta_dy"] * mesh_data["dgam_dx"]
        )
    )

    # del mesh_data["dxsi_dx"]
    # del mesh_data["dxsi_dy"]
    # del mesh_data["dxsi_dz"]
    # del mesh_data["deta_dx"]
    # del mesh_data["deta_dy"]
    # del mesh_data["deta_dz"]
    # del mesh_data["dgam_dx"]
    # del mesh_data["dgam_dy"]
    # del mesh_data["dgam_dz"]

    # Fortran array index starts at 1, now convert to C array index beginning at 0
    mesh_data["ibool"] = mesh_data["ibool"] - 1

    # use xyz_glob
    nglob = mesh_data["nglob"]
    xyz_glob = np.zeros((nglob, 3))
    xyz_glob[:, 0] = mesh_data["x"]
    xyz_glob[:, 1] = mesh_data["y"]
    xyz_glob[:, 2] = mesh_data["z"]
    mesh_data["xyz_glob"] = xyz_glob

    del mesh_data["x"]
    del mesh_data["y"]
    del mesh_data["z"]

    # add xyz_elem
    iglob_elem = mesh_data["ibool"][:, MIDZ, MIDY, MIDX]
    mesh_data["xyz_elem"] = mesh_data["xyz_glob"][iglob_elem, :]

    # FIXME bad idea due to 410 undulation. Need to modify the specfem code
    ## separate mesh layers across 410-km
    ## 40: above 410, 41: below 410
    # idoubling = mesh_data['idoubling']
    # depth = (1.0 - np.sum(mesh_data['xyz_elem']**2, axis=0)*0.5) * R_EARTH_KM

    ## this is dangerous due to 410 undulation
    # ii = (idoubling == IFLAG_670_220) & (depth < 410)
    # idoubling[ii] = 10*IFLAG_670_220

    # ii = (idoubling == IFLAG_670_220) & (depth > 410)
    # idoubling[ii] = 10*IFLAG_670_220 + 1

    # nspec = int(mesh_data['nspec'])
    #  for ispec in range(nspec):
    #    for i in range(NGLLX):
    #      for j in range(NGLLY):
    #        for k in range(NGLLZ):
    #          iglob = mesh_data['ibool'][i,j,k,ispec] - 1
    #          xyz_gll[0,i,j,k,ispec] = mesh_data['x'][iglob]
    #          xyz_gll[1,i,j,k,ispec] = mesh_data['y'][iglob]
    #          xyz_gll[2,i,j,k,ispec] = mesh_data['z'][iglob]
    #  xyz_gll = np.zeros((3,NGLLX,NGLLY,NGLLZ,nspec))

    return mesh_data


def sem_mesh_mpi_read(mesh_mpi_file):
    """read in SEM solver_data_mpi.bin"""
    from scipy.io import FortranFile

    mesh_mpi_data = {}

    with FortranFile(mesh_mpi_file, "r") as f:
        for field in MESH_MPI_ARRAY_LIST:
            field_name = field[0]
            data_type = field[1]
            mesh_mpi_data[field_name] = f.read_ints(dtype=data_type)

    num_interfaces = mesh_mpi_data["num_interfaces"][0]
    mesh_mpi_data["num_interfaces"] = num_interfaces

    max_nibool_interfaces = mesh_mpi_data["max_nibool_interfaces"][0]
    mesh_mpi_data["max_nibool_interfaces"] = max_nibool_interfaces

    dims = (
        num_interfaces,
        max_nibool_interfaces,
    )  # inversed order from Fortran convention
    mesh_mpi_data["ibool_interfaces"] = np.reshape(
        mesh_mpi_data["ibool_interfaces"] - 1, dims
    )

    return mesh_mpi_data


def sem_mesh_get_vol_gll(mesh_data):
    """get xyz and volume weights of each gll point"""

    from gll_library import zwgljd

    # --- quadrature weights on GLL points
    zx, wx = zwgljd(NGLLX, GAUSSALPHA, GAUSSBETA)
    zy, wy = zwgljd(NGLLY, GAUSSALPHA, GAUSSBETA)
    zz, wz = zwgljd(NGLLZ, GAUSSALPHA, GAUSSBETA)

    # wgll_cube = wx.reshape((NGLLX,1,1))*wy.reshape((1,NGLLY,1))*wx.reshape((1,1,NGLLZ))

    # --- jacobian * gll_quad_weights
    # vol_gll = mesh_data['jacobian']*wgll_cube.reshape((NGLLX,NGLLY,NGLLZ,1))
    # vol_gll = np.array(mesh_data['jacobian']*wx[:,None,None,None]*wy[None,:,None,None]*wz[None,None,:,None], dtype='float32')
    vol_gll = (
        mesh_data["jacobian"]
        * wz[None, :, None, None]
        * wy[None, None, :, None]
        * wx[None, None, None, :]
    )

    return vol_gll


# @numba.jit("void(float64[:,:], float64[:], float64[:], float64[:,:])")
# def sem_jacobian_hex27(anchors_xyz, uvw, xyz, dudx):
# @numba.jit
# @numba.jit("void(float64[:,::1], float64[::1], float64[::1], float64[:,::1])")
@numba.jit(nogil=True) #, cache=True)
def sem_jacobian_hex27(anchors_xyz, uvw, xyz, dudx):
    """
    compute 3D jacobian at a given point for a 27-node element
    map from local coordinate (uvw) to physical position (xyz)
    the shape the element is defined by the anchor points (xyz_anchor)
    !
    !-input
    ! anchors_xyz(27,3): xyz of anchor points, the order of the 27 nodes must be from
    !   subroutine anchor_index_hex27()
    ! uvw(3): local coordinate
    !
    !-output
    ! xyz(3): map uvw to physical space
    ! dudx(3,3): jacobian matrix
    """
    # lagrange basis at nodes [-1,0,1] evaluated at uvw
    lag = np.zeros((3, 3))
    lag[0, :] = uvw * (uvw - 1.0) / 2.0
    lag[1, :] = 1.0 - uvw**2
    lag[2, :] = uvw * (uvw + 1.0) / 2.0

    # derivative of lagrange polynomials at uvw
    dlag = np.zeros((3, 3))
    dlag[0, :] = uvw - 0.5
    dlag[1, :] = -2.0 * uvw
    dlag[2, :] = uvw + 0.5

    # HEX27_IJK = np.array(
    #     [
    #         # 8 corners
    #         [0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0], [0, 0, 2], [2, 0, 2], [2, 2, 2], [0, 2, 2],
    #         # 16 edge centers
    #         [1, 0, 0], [2, 1, 0], [1, 2, 0], [0, 1, 0], [0, 0, 1], [2, 0, 1], [2, 2, 1], [0, 2, 1],
    #         [1, 0, 2], [2, 1, 2], [1, 2, 2], [0, 1, 2],
    #         # 6 face centers
    #         [1, 1, 0], [1, 0, 1], [2, 1, 1], [1, 2, 1], [0, 1, 1], [1, 1, 2],
    #         # 1 body center
    #         [1, 1, 1],
    #     ],
    #     dtype=np.int64,
    # )
    ii, jj, kk = ANCHOR_NODES[:, 0], ANCHOR_NODES[:, 1], ANCHOR_NODES[:, 2]
    shape3D = lag[ii, 0] * lag[jj, 1] * lag[kk, 2]
    dershape3D = np.zeros((3, 27))
    dershape3D[0, :] = dlag[ii, 0] *  lag[jj, 1] *  lag[kk, 2]
    dershape3D[1, :] =  lag[ii, 0] * dlag[jj, 1] *  lag[kk, 2]
    dershape3D[2, :] =  lag[ii, 0] *  lag[jj, 1] * dlag[kk, 2]

    # for a in range(27):
    #     i, j, k = ANCHOR_NODES[a, :]
    #     # xyz += anchors_xyz[a,:] * (lag[i, 0] * lag[j, 1] * lag[k, 2])
    #     # dxdu[0, :] += anchors_xyz[a, 0] * (dlag[i,0] * lag[j,1] * lag[k, 2])
    #     # dxdu[1, :] += anchors_xyz[a, 1] * (lag[i,0] * dlag[j,1] * lag[k, 2])
    #     # dxdu[2, :] += anchors_xyz[a, 2] * (lag[i,0] * lag[j,1] * dlag[k, 2])
    #     shape3D[a] = lag[i, 0] * lag[j, 1] * lag[k, 2]
    #     dershape3D[a, 0] = dlag[i, 0] * lag[j, 1] * lag[k, 2]
    #     dershape3D[a, 1] = lag[i, 0] * dlag[j, 1] * lag[k, 2]
    #     dershape3D[a, 2] = lag[i, 0] * lag[j, 1] * dlag[k, 2]

    # xyz[:] = 0
    # xyz = np.zeros(3)
    # dxdu = np.zeros((3, 3))
    xyz[:] = np.dot(shape3D, anchors_xyz)
    dxdu = np.transpose(np.dot(dershape3D, anchors_xyz))
    # for a in range(27):
    #     xyz += shape3D[a] * anchors_xyz[a, :]
    #     dxdu[0, :] += anchors_xyz[a, 0] * dershape3D[a, :]
    #     dxdu[1, :] += anchors_xyz[a, 1] * dershape3D[a, :]
    #     dxdu[2, :] += anchors_xyz[a, 2] * dershape3D[a, :]

    # dudx = np.linalg.inv(dxdu)

    # transpose(adjoint(dxdu))
    # dudx = np.zeros((3, 3))
    dudx[0, 0] = dxdu[1, 1] * dxdu[2, 2] - dxdu[2, 1] * dxdu[1, 2]
    dudx[1, 0] = -(dxdu[1, 0] * dxdu[2, 2] - dxdu[2, 0] * dxdu[1, 2])
    dudx[2, 0] = dxdu[1, 0] * dxdu[2, 1] - dxdu[2, 0] * dxdu[1, 1]
    dudx[0, 1] = -(dxdu[0, 1] * dxdu[2, 2] - dxdu[2, 1] * dxdu[0, 2])
    dudx[1, 1] = dxdu[0, 0] * dxdu[2, 2] - dxdu[2, 0] * dxdu[0, 2]
    dudx[2, 1] = -(dxdu[0, 0] * dxdu[2, 1] - dxdu[2, 0] * dxdu[0, 1])
    dudx[0, 2] = dxdu[0, 1] * dxdu[1, 2] - dxdu[1, 1] * dxdu[0, 2]
    dudx[1, 2] = -(dxdu[0, 0] * dxdu[1, 2] - dxdu[1, 0] * dxdu[0, 2])
    dudx[2, 2] = dxdu[0, 0] * dxdu[1, 1] - dxdu[1, 0] * dxdu[0, 1]
    # det(dxdu)
    det = dxdu[0, 0] * dudx[0, 0] + dxdu[0, 1] * dudx[1, 0] + dxdu[0, 2] * dudx[2, 0]
    # inverse of dxdu: transpose(adjoint(dxdu)) / det(dxdu)
    dudx[:] = dudx[:] / det

    # return xyz, dudx


# @numba.njit("Tuple((float64, boolean))(float64[:,:], float64[:], float64[:], int64)")
# @numba.jit()
# @numba.njit("Tuple((float64, boolean))(float64[:,::1], float64[::1], float64[::1], int64)")
@numba.jit(nogil=True) #, cache=True)
def sem_map2cube_hex27(anchors_xyz, target_xyz, located_uvw, max_niter=5):
    """
    !-map a given point in physical space (xyz) to the
    ! reference cube (uvw) for a 27-node hexahedron element (8 vertex, 12 egde
    ! centers, 6 face centers and 1 body center),
    ! and also flag whether the point is inside the cube
    ! if the point lies outside the element, calculate the bounded (xi,eta,gamma)
    ! inside or on the surface of the reference unit cube.
    !
    !-inputs:
    ! (real) xyz_anchor(27,3): anchor points of the element
    ! (real) xyz(3): coordinates of the target point
    !
    !-outputs:
    ! (real) uvw(3): local coordinates in reference cube
    ! (real) misloc: location misfit abs(xyz - XYZ(uvw))
    ! (logical) flag_inside: flag whether the target point locates inside the element
    """
    # find nearest anchor nodes as initial uvw
    # min_dist_sq = 1.0e5
    # best_a = 0
    # for a in range(27):
    #     dist_sq = sum((anchors_xyz[a, :] - target_xyz[:]) ** 2)
    #     if dist_sq < min_dist_sq:
    #         min_dist_sq = dist_sq
    #         best_a = a
    # i, j, k = ANCHOR_GLL_INDEX[best_a, :]
    # located_uvw[0] = GLL_NODES[i]
    # located_uvw[1] = GLL_NODES[j]
    # located_uvw[2] = GLL_NODES[k]

    # TODO the above found initial value actually increases final mis-location???
    located_uvw[:] = 0.0

    # iteratively update local coordinate uvw to approach the target xyz
    xyz = np.zeros(3)
    dudx = np.zeros((3,3))
    is_inside = True
    for iter in range(max_niter):
        # predicted xyzi and Jacobian for the current uvw
        sem_jacobian_hex27(anchors_xyz, located_uvw, xyz, dudx)
        # xyz, dudx = sem_jacobian_hex27(anchors_xyz, located_uvw)
        # print(f"{iter=}, {located_uvw=}, {xyz=}")

        # update
        located_uvw[:] = located_uvw[:] + np.dot(dudx, target_xyz - xyz)

        is_inside = True
        mask = located_uvw < -1
        if np.any(mask): is_inside = False
        located_uvw[mask] = -1
        mask = located_uvw > 1
        if np.any(mask): is_inside = False
        located_uvw[mask] = 1

        # # limit inside the reference cube
        # is_inside = True
        # for i in range(3):
        #     if located_uvw[i] < -1:
        #         located_uvw[i] = -1
        #         is_inside = False
        #     elif located_uvw[i] > 1:
        #         located_uvw[i] = 1
        #         is_inside = False

    sem_jacobian_hex27(anchors_xyz, located_uvw, xyz, dudx)
    # xyz, dudx = sem_jacobian_hex27(anchors_xyz, located_uvw)
    misloc = (sum(target_xyz - xyz) ** 2) ** 0.5
    # print(f"{iter=}, {located_uvw=}, {xyz=}")

    return misloc, is_inside


def sem_mesh_locate_points(
    mesh_data, xyz, idoubling=-1, kdtree_num_element=2.0, max_dist_ratio=2.0
):
    """locate points in the SEM mesh.
    mesh_data: return value from sem_mesh_read()
    xyz(n,3): xyz of n points
    idoubling: integer or integer array of size n. idoubling for the n points which denotes mesh regions (surface-Moho,Moho-410,410-660,etc)
      -1 means no specific region and interpolation will be done to all the elements in the mesh, otherwise interpolation is only done for those elements with the same idoubling value.
    kdtree_num_element: radius factor as number of multiples of the maximum element half size used in kdtree search of neighboring elements to target point.
    max_dist_ratio: maximum ratio between the distance from target point to the element center and the element half size. Used to ignore element which is too far away from the target point. Sometimes if this value is too close to one, the target point slightly outside the mesh will be marked as NOT located, even if the SEM could allow a point outside the element be located.

    output:
      status(n): -1=not located,0=close to but outside the element,1=inside element
      ispec(n): element num that located
      uvw(n,3): local coordinates of n points
      misloc(n): location residual
      misratio(n): misloc/element_half_size
    """
    from scipy import spatial

    # from gll_library import zwgljd, lagrange_poly
    # from jacobian_hex27 import xyz2cube_bounded_hex27, anchor_index_hex27

    if max_dist_ratio < 1:
        warnings.warn(
            "max_dist_ratio should be larger than one! Default value 2.0 will be used."
        )
        max_dist_ratio = 2.0

    npoints = xyz.shape[0]
    idoubling = np.array(idoubling, dtype="int")
    if idoubling.size == 1:
        idoubling = np.ones(npoints, dtype="int") * int(idoubling)
    elif idoubling.size != npoints:
        raise Exception(
            "idoubling must either be an integer or an integer array of npoints "
        )

    nspec = mesh_data["nspec"]
    ibool = mesh_data["ibool"]
    source_idoubling = mesh_data["idoubling"]
    xyz_glob = mesh_data["xyz_glob"]
    xyz_elem = mesh_data["xyz_elem"]

    # --- kdtree search nearby elements around each target point
    tree_elem = spatial.cKDTree(
        np.column_stack((xyz_elem[:, 0], xyz_elem[:, 1], xyz_elem[:, 2]))
    )

    tree_xyz = spatial.cKDTree(np.column_stack((xyz[:, 0], xyz[:, 1], xyz[:, 2])))

    # determine element size (approximately)
    element_half_size = np.zeros(nspec)
    for ispec in range(nspec):
        # distance between gll points and the central gll point
        iglob1 = ibool[ispec, :, :, :].ravel()
        dist = (
            np.sum((xyz_elem[ispec, :] - xyz_glob[iglob1, :]) ** 2, axis=1)
            ** 0.5
        )
        element_half_size[ispec] = np.max(dist)

    # get neighbouring elements around each target location xyz
    neighbor_lists = tree_xyz.query_ball_tree(
        tree_elem, kdtree_num_element * np.max(element_half_size)
    )

    # --- loop over each point, get the location info
    iax = ANCHOR_GLL_INDEX[:,0]
    iay = ANCHOR_GLL_INDEX[:,1]
    iaz = ANCHOR_GLL_INDEX[:,2]

    # xigll, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
    # yigll, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
    # zigll, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

    status_all = np.zeros(npoints, dtype="int")
    status_all[:] = -1
    ispec_all = np.zeros(npoints, dtype="int")
    uvw_all = np.zeros((npoints, 3))
    misloc_all = np.zeros(npoints)
    misloc_all[:] = np.inf
    misratio_all = np.zeros(npoints)

    ipoint_select = [ipoint for ipoint in range(npoints) if neighbor_lists[ipoint]]
    uvw = np.zeros(3)
    # for ipoint in range(npoints):
    for ipoint in ipoint_select:
        # if not neighbor_lists[ipoint]: continue
        # get neibouring elements
        # convert list to numpy array to have index slicing
        ispec_list = np.array(neighbor_lists[ipoint])
        # ratio between the distance from target point to the element center and the element size
        # print(ispec_list.shape)
        # print(xyz_elem.shape)
        # print(xyz.shape)
        dist_ratio = (
            np.sum((xyz_elem[ispec_list, :] - xyz[ipoint, :]) ** 2, axis=1)
            ** 0.5
            / element_half_size[ispec_list]
        )
        # remove elements too far away from target point
        idx = dist_ratio < max_dist_ratio
        # skip elements that does NOT have the same idoubling as xyz
        if idoubling[ipoint] != -1:
            idx = idx & (source_idoubling[ispec_list] == idoubling[ipoint])
        if not np.any(idx):
            continue
        ispec_list = ispec_list[idx]
        dist_ratio = dist_ratio[idx]
        # loop each element, start from the closest element
        for ispec in ispec_list[np.argsort(dist_ratio)]:
            # if (idoubling[ipoint] != -1 and
            #    idoubling[ipoint] != source_idoubling[ispec]):
            #  continue
            iglob = ibool[ispec, iaz, iay, iax]
            xyz_anchor = xyz_glob[iglob, :]
            misloc, is_inside = sem_map2cube_hex27(xyz_anchor, xyz[ipoint, :], uvw)
            ##DEBUG
            # if is_inside and status_all[ipoint]==1:
            #  warnings.warn("point is located inside more than one element",
            #      xyz[:,ipoint], xyz_anchor)
            if misloc > misloc_all[ipoint] and is_inside:
                warnings.warn(
                    "point located inside an element but with a larger misloc: current/previous = %f/%f"
                    % (misloc, misloc_all[ipoint])
                )
            if misloc < misloc_all[ipoint] or is_inside:
                status_all[ipoint] = 0
                ispec_all[ipoint] = ispec
                uvw_all[ipoint, :] = uvw
                misloc_all[ipoint] = misloc
                misratio_all[ipoint] = misloc / element_half_size[ispec]
            # skip the rest elements since points already located inside an element
            # this means if multiple elements overlap (should not occur) we only take the first found element where the point locates inside
            if is_inside:
                status_all[ipoint] = 1
                break
        # if 'uvw' in loc_data[ipoint]:
        #  hlagx = lagrange_poly(xigll, uvw[0])
        #  hlagy = lagrange_poly(yigll, uvw[1])
        #  hlagz = lagrange_poly(zigll, uvw[2])
        #  loc_data[ipoint]['lagrange'] = hlagx[:,None,None]*hlagy[None,:,None]*hlagz[None,None,:]

    return status_all, ispec_all, uvw_all, misloc_all, misratio_all


def sem_boundary_mesh_read(mesh_file):
    """read in SEM mesh slice"""
    from scipy.io import FortranFile

    mesh_data = {}

    with FortranFile(mesh_file, "r") as f:
        for field in BOUNDARY_ARRAY_LIST:
            field_name = field[0]
            data_type = field[1]
            mesh_data[field_name] = f.read_ints(dtype=data_type)

    mesh_data["nspec2D_teleseismic_xmin"] = mesh_data["nspec2D_teleseismic_xmin"][0]
    mesh_data["nspec2D_teleseismic_xmax"] = mesh_data["nspec2D_teleseismic_xmax"][0]
    mesh_data["nspec2D_teleseismic_ymin"] = mesh_data["nspec2D_teleseismic_ymin"][0]
    mesh_data["nspec2D_teleseismic_ymax"] = mesh_data["nspec2D_teleseismic_ymax"][0]
    mesh_data["nspec2D_teleseismic_zmin"] = mesh_data["nspec2D_teleseismic_zmin"][0]

    # reshape
    # NB: binary files are written in Fortran column-major convention !!!
    # NB: reshape 1-D array to matrix by Fortran convention, and
    # NB: also convert to a contiguous array in memory, in case of direct memory copy or MPI transfer
    mesh_data["area_teleseismic_xmin"] = np.ascontiguousarray(
        np.reshape(mesh_data["area_teleseismic_xmin"], (NGLLY, NGLLZ, -1), order="F")
    )
    mesh_data["area_teleseismic_xmax"] = np.ascontiguousarray(
        np.reshape(mesh_data["area_teleseismic_xmax"], (NGLLY, NGLLZ, -1), order="F")
    )
    mesh_data["area_teleseismic_ymin"] = np.ascontiguousarray(
        np.reshape(mesh_data["area_teleseismic_ymin"], (NGLLX, NGLLZ, -1), order="F")
    )
    mesh_data["area_teleseismic_ymax"] = np.ascontiguousarray(
        np.reshape(mesh_data["area_teleseismic_ymax"], (NGLLX, NGLLZ, -1), order="F")
    )
    mesh_data["area_teleseismic_zmin"] = np.ascontiguousarray(
        np.reshape(mesh_data["area_teleseismic_zmin"], (NGLLX, NGLLY, -1), order="F")
    )

    # cut data arrays to lengths actually used
    mesh_data["ibelm_teleseismic_xmin"] = mesh_data["ibelm_teleseismic_xmin"][
        0 : mesh_data["nspec2D_teleseismic_xmin"]
    ]
    mesh_data["ibelm_teleseismic_xmax"] = mesh_data["ibelm_teleseismic_xmax"][
        0 : mesh_data["nspec2D_teleseismic_xmax"]
    ]
    mesh_data["ibelm_teleseismic_ymin"] = mesh_data["ibelm_teleseismic_ymin"][
        0 : mesh_data["nspec2D_teleseismic_ymin"]
    ]
    mesh_data["ibelm_teleseismic_ymax"] = mesh_data["ibelm_teleseismic_ymax"][
        0 : mesh_data["nspec2D_teleseismic_ymax"]
    ]
    mesh_data["ibelm_teleseismic_zmin"] = mesh_data["ibelm_teleseismic_zmin"][
        0 : mesh_data["nspec2D_teleseismic_zmin"]
    ]

    mesh_data["area_teleseismic_xmin"] = mesh_data["area_teleseismic_xmin"][
        :, :, 0 : mesh_data["nspec2D_teleseismic_xmin"]
    ]
    mesh_data["area_teleseismic_xmax"] = mesh_data["area_teleseismic_xmax"][
        :, :, 0 : mesh_data["nspec2D_teleseismic_xmax"]
    ]
    mesh_data["area_teleseismic_ymin"] = mesh_data["area_teleseismic_ymin"][
        :, :, 0 : mesh_data["nspec2D_teleseismic_ymin"]
    ]
    mesh_data["area_teleseismic_ymax"] = mesh_data["area_teleseismic_ymax"][
        :, :, 0 : mesh_data["nspec2D_teleseismic_ymax"]
    ]
    mesh_data["area_teleseismic_zmin"] = mesh_data["area_teleseismic_zmin"][
        :, :, 0 : mesh_data["nspec2D_teleseismic_zmin"]
    ]

    return mesh_data


# @numba.jit
def assemble_MPI_scalar(
    array_glob,
    num_interfaces,
    max_nibool_interfaces,
    nibool_interfaces,
    ibool_interfaces,
    my_neighbors,
):
    from mpi4py import MPI

    # import numba, numba_mpi
    # integer, intent(in) :: num_interfaces,max_nibool_interfaces
    # integer, dimension(num_interfaces), intent(in) :: nibool_interfaces,my_neighbors
    # integer, dimension(max_nibool_interfaces,num_interfaces), intent(in) :: ibool_interfaces

    comm = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()
    # mpi_rank = numba_mpi.rank()

    buffer_send_scalar = np.zeros(
        (num_interfaces, max_nibool_interfaces), dtype=array_glob.dtype
    )
    buffer_recv_scalar = np.zeros(
        (num_interfaces, max_nibool_interfaces), dtype=array_glob.dtype
    )

    for iinterface in range(num_interfaces):
        npts = nibool_interfaces[iinterface]
        buffer_send_scalar[iinterface, :npts] = array_glob[
            ibool_interfaces[iinterface, :npts]
        ]

    # send messages
    # send_requests = np.zeros((num_interfaces,), dtype=numba_mpi.RequestType)
    # recv_requests = np.zeros((num_interfaces,), dtype=numba_mpi.RequestType)
    send_requests = []
    recv_requests = []
    for iinterface in range(num_interfaces):
        npts = nibool_interfaces[iinterface]
        req = comm.Isend(
            # status, req = numba_mpi.isend(
            buffer_send_scalar[iinterface, :npts],
            dest=my_neighbors[iinterface],
            tag=11,
        )
        send_requests.append(req)
        # send_requests[iinterface] = req[0]
        req = comm.Irecv(
            # status, req = numba_mpi.irecv(
            buffer_recv_scalar[iinterface, :npts],
            source=my_neighbors[iinterface],
            tag=11,
        )
        recv_requests.append(req)
        # recv_requests[iinterface] = req[0]

    # wait for communications completion (recv)
    MPI.Request.Waitall(recv_requests)
    # numba_mpi.waitall(recv_requests)

    # adding contributions of neighbors, in the order of mpi_ranks
    sum_glob = np.zeros_like(array_glob)

    inds = np.argsort(my_neighbors)
    mask = my_neighbors[inds] < mpi_rank

    for i in inds[mask]:
        npts = nibool_interfaces[i]
        sum_glob[ibool_interfaces[i, :npts]] += buffer_recv_scalar[i, :npts]

    sum_glob[:] += array_glob[:]

    for i in inds[~mask]:
        npts = nibool_interfaces[i]
        sum_glob[ibool_interfaces[i, :npts]] += buffer_recv_scalar[i, :npts]

    # for i in range(num_interfaces):
    #     npts = nibool_interfaces[i]
    #     array_glob[ibool_interfaces[i, :npts]] += buffer_recv_scalar[i, :npts]

    # wait for communications completion (send)
    MPI.Request.Waitall(send_requests)
    # numba_mpi.waitall(send_requests)

    array_glob[:] = sum_glob[:]


@numba.jit(nogil=True) #, cache=True)
def gll2glob(
    u_gll,
    nglob,
    ibool,  # [0:nspec,0:nzgll,0:nygll,0:nxgll]
):
    """
    u_gll[nspec,nzgll,nygll,nxgll]
    """
    nspec, ngllz, nglly, ngllx = ibool.shape
    u_glob = np.zeros(nglob, dtype=u_gll.dtype)
    # counts = np.zeros(nglob, dtype=ibool.dtype)

    for e in range(nspec):
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    # counts[idof] += 1
                    u_glob[idof] = u_glob[idof] + u_gll[e, k, j, i]
    # u_glob = u_glob / counts
    return u_glob


@numba.jit(nogil=True) #, cache=True)
def laplacian_iso(
    u_glob,
    kappa,
    wgll,
    dlag_gll,
    ibool,
    dxsi_dx,
    dxsi_dy,
    dxsi_dz,
    deta_dx,
    deta_dy,
    deta_dz,
    dgam_dx,
    dgam_dy,
    dgam_dz,
    jacobian,
):
    """
    calculate weak form <phi_gll, grad(K * grad(u_glob))> = -1 * int(grad(phi_gll) * K * grad(u_glob), dV)
    u_glob -> u_{npsec,ngllz,nglly,ngllx} is the trial function
    phi_gll_{npsec,ngllz,nglly,ngllx} is the test function

    NOTE: need to assemble the results if run in parralel (assemble_MPI_scalar)
    """
    nspec, ngllz, nglly, ngllx = ibool.shape
    # assert u_glob.shape == out_glob.shape

    dtype = u_glob.dtype
    sl = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_xsil = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_etal = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_gaml = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    stif = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    out_glob = np.zeros(u_glob.shape, dtype=dtype)

    for e in range(nspec):

        # sl = u_glob[ibool[e, :, :, :]]
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    sl[k, j, i] = u_glob[idof]

        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Compute derivatives at elemental level
                    dsl_dxsi = 0
                    dsl_deta = 0
                    dsl_dgam = 0

                    for l in range(ngllx):
                        dsl_dxsi = dsl_dxsi + sl[k, j, l] * dlag_gll[l, i]
                        dsl_deta = dsl_deta + sl[k, l, i] * dlag_gll[l, j]
                        dsl_dgam = dsl_dgam + sl[l, j, i] * dlag_gll[l, k]

                    # Get derivatives informations
                    dxsi_dxl = dxsi_dx[e, k, j, i]
                    dxsi_dyl = dxsi_dy[e, k, j, i]
                    dxsi_dzl = dxsi_dz[e, k, j, i]
                    deta_dxl = deta_dx[e, k, j, i]
                    deta_dyl = deta_dy[e, k, j, i]
                    deta_dzl = deta_dz[e, k, j, i]
                    dgam_dxl = dgam_dx[e, k, j, i]
                    dgam_dyl = dgam_dy[e, k, j, i]
                    dgam_dzl = dgam_dz[e, k, j, i]
                    jacobianl = jacobian[e, k, j, i]

                    # Get physical derivatives
                    dsl_dxl = (
                        dsl_dxsi * dxsi_dxl + dsl_deta * deta_dxl + dsl_dgam * dgam_dxl
                    )
                    dsl_dyl = (
                        dsl_dxsi * dxsi_dyl + dsl_deta * deta_dyl + dsl_dgam * dgam_dyl
                    )
                    dsl_dzl = (
                        dsl_dxsi * dxsi_dzl + dsl_deta * deta_dzl + dsl_dgam * dgam_dzl
                    )

                    # Start product: Kij * du_dxi * dxsi_dxj
                    # kl = kappa_gll[e, k, j, i]
                    grad_xsil[k, j, i] = (
                        jacobianl
                        * kappa
                        * (dxsi_dxl * dsl_dxl + dxsi_dyl * dsl_dyl + dxsi_dzl * dsl_dzl)
                    )
                    grad_etal[k, j, i] = (
                        jacobianl
                        * kappa
                        * (deta_dxl * dsl_dxl + deta_dyl * dsl_dyl + deta_dzl * dsl_dzl)
                    )
                    grad_gaml[k, j, i] = (
                        jacobianl
                        * kappa
                        * (dgam_dxl * dsl_dxl + dgam_dyl * dsl_dyl + dgam_dzl * dsl_dzl)
                    )

        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Compute derivatives at elemental level
                    lapla_x = 0.0
                    lapla_y = 0.0
                    lapla_z = 0.0
                    for l in range(ngllx):
                        lapla_x = (
                            lapla_x + grad_xsil[k, j, l] * dlag_gll[i, l] * wgll[l]
                        )
                        lapla_y = (
                            lapla_y + grad_etal[k, l, i] * dlag_gll[j, l] * wgll[l]
                        )
                        lapla_z = (
                            lapla_z + grad_gaml[l, j, i] * dlag_gll[k, l] * wgll[l]
                        )

                    # Stiffness
                    stif[k, j, i] = (
                        wgll[j] * wgll[k] * lapla_x
                        + wgll[i] * wgll[k] * lapla_y
                        + wgll[i] * wgll[j] * lapla_z
                    )

        # Go back to dof
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    out_glob[idof] += stif[k, j, i]

    return -1 * out_glob


@numba.jit(nogil=True) #, cache=True)
def laplacian_iso3D(
    u_glob,
    kappa_gll,
    wgll,
    dlag_gll,
    ibool,
    dxsi_dx,
    dxsi_dy,
    dxsi_dz,
    deta_dx,
    deta_dy,
    deta_dz,
    dgam_dx,
    dgam_dy,
    dgam_dz,
    jacobian,
):
    """
    int(grad(phi_gll) * K * grad(u), dV)
    u_glob -> u_{npsec,ngllz,nglly,ngllx} is the trial function
    phi_gll_{npsec,ngllz,nglly,ngllx} is the test function

    NOTE: need to assemble the results if run in parralel (assemble_MPI_scalar)
    """
    nspec, ngllz, nglly, ngllx = ibool.shape
    # assert u_glob.shape == out_glob.shape

    dtype = u_glob.dtype
    sl = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_xsil = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_etal = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    grad_gaml = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    stif = np.zeros((ngllz, nglly, ngllx), dtype=dtype)
    out_glob = np.zeros(u_glob.shape, dtype=dtype)

    for e in range(nspec):

        # sl = u_glob[ibool[e, :, :, :]]
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    sl[k, j, i] = u_glob[idof]

        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Compute derivatives at elemental level
                    dsl_dxsi = 0
                    dsl_deta = 0
                    dsl_dgam = 0

                    for l in range(ngllx):
                        dsl_dxsi = dsl_dxsi + sl[k, j, l] * dlag_gll[l, i]
                        dsl_deta = dsl_deta + sl[k, l, i] * dlag_gll[l, j]
                        dsl_dgam = dsl_dgam + sl[l, j, i] * dlag_gll[l, k]

                    # Get derivatives informations
                    dxsi_dxl = dxsi_dx[e, k, j, i]
                    dxsi_dyl = dxsi_dy[e, k, j, i]
                    dxsi_dzl = dxsi_dz[e, k, j, i]
                    deta_dxl = deta_dx[e, k, j, i]
                    deta_dyl = deta_dy[e, k, j, i]
                    deta_dzl = deta_dz[e, k, j, i]
                    dgam_dxl = dgam_dx[e, k, j, i]
                    dgam_dyl = dgam_dy[e, k, j, i]
                    dgam_dzl = dgam_dz[e, k, j, i]
                    jacobianl = jacobian[e, k, j, i]

                    # Get physical derivatives
                    dsl_dxl = (
                        dsl_dxsi * dxsi_dxl + dsl_deta * deta_dxl + dsl_dgam * dgam_dxl
                    )
                    dsl_dyl = (
                        dsl_dxsi * dxsi_dyl + dsl_deta * deta_dyl + dsl_dgam * dgam_dyl
                    )
                    dsl_dzl = (
                        dsl_dxsi * dxsi_dzl + dsl_deta * deta_dzl + dsl_dgam * dgam_dzl
                    )

                    # Start product: Kij * du_dxi * dxsi_dxj
                    kappal = kappa_gll[e, k, j, i]
                    grad_xsil[k, j, i] = (
                        jacobianl
                        * kappal
                        * (dxsi_dxl * dsl_dxl + dxsi_dyl * dsl_dyl + dxsi_dzl * dsl_dzl)
                    )
                    grad_etal[k, j, i] = (
                        jacobianl
                        * kappal
                        * (deta_dxl * dsl_dxl + deta_dyl * dsl_dyl + deta_dzl * dsl_dzl)
                    )
                    grad_gaml[k, j, i] = (
                        jacobianl
                        * kappal
                        * (dgam_dxl * dsl_dxl + dgam_dyl * dsl_dyl + dgam_dzl * dsl_dzl)
                    )

        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Compute derivatives at elemental level
                    lapla_x = 0.0
                    lapla_y = 0.0
                    lapla_z = 0.0
                    for l in range(ngllx):
                        lapla_x = (
                            lapla_x + grad_xsil[k, j, l] * dlag_gll[i, l] * wgll[l]
                        )
                        lapla_y = (
                            lapla_y + grad_etal[k, l, i] * dlag_gll[j, l] * wgll[l]
                        )
                        lapla_z = (
                            lapla_z + grad_gaml[l, j, i] * dlag_gll[k, l] * wgll[l]
                        )

                    # Stiffness
                    stif[k, j, i] = (
                        wgll[j] * wgll[k] * lapla_x
                        + wgll[i] * wgll[k] * lapla_y
                        + wgll[i] * wgll[j] * lapla_z
                    )

        # Go back to dof
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    idof = ibool[e, k, j, i]
                    out_glob[idof] += stif[k, j, i]

    return -1 * out_glob
