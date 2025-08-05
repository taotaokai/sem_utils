#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

"""
parameters from constants.h
"""

# -- radius of the Earth
R_EARTH = 6371000.0
R_EARTH_KM = 6371.0

# --- GLL

#! number of GLL points in each direction of an element (degree plus one)
NGLLX = 5
NGLLY = NGLLX
NGLLZ = NGLLX

GLL5_NODES = np.array([-1, -((3.0 / 7.0) ** 0.5), 0, (3.0 / 7.0) ** 0.5, 1])
GLL5_WEIGHTS = np.array([0.1, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 0.1])

#! mid-points inside a GLL element
MIDX = NGLLX // 2
MIDY = NGLLY // 2
MIDZ = NGLLZ // 2

#! for the Gauss-Lobatto-Legendre points and weights
GAUSSALPHA = 0
GAUSSBETA = 0

# --- define flag for elements
IFLAG_CRUST = 1

IFLAG_80_MOHO = 2
IFLAG_220_80 = 3
IFLAG_670_220 = 4
IFLAG_MANTLE_NORMAL = 5

IFLAG_OUTER_CORE_NORMAL = 6

IFLAG_INNER_CORE_NORMAL = 7
IFLAG_MIDDLE_CENTRAL_CUBE = 8
IFLAG_BOTTOM_CENTRAL_CUBE = 9
IFLAG_TOP_CENTRAL_CUBE = 10
IFLAG_IN_FICTITIOUS_CUBE = 11

IFLAG_DUMMY = 100

# --- data structure of proc000*_reg1_solver_data.bin
MESH_ARRAY_LIST = [
    ("nspec", "i4"),
    ("nglob", "i4"),
    ("x", "f4"),
    ("y", "f4"),
    ("z", "f4"),
    ("ibool", "i4"),
    ("idoubling", "i4"),
    ("ispec_is_tiso", "i4"),
    ("dxsi_dx", "f4"),
    ("dxsi_dy", "f4"),
    ("dxsi_dz", "f4"),
    ("deta_dx", "f4"),
    ("deta_dy", "f4"),
    ("deta_dz", "f4"),
    ("dgam_dx", "f4"),
    ("dgam_dy", "f4"),
    ("dgam_dz", "f4"),
]

# --- data structure of proc000*_reg1_solver_data.bin
# src/meshfem3D/save_arrays_solver.f90: subroutine save_MPI_arrays
#  integer, intent(in) :: num_interfaces,max_nibool_interfaces
#  integer, dimension(num_interfaces), intent(in) :: nibool_interfaces,my_neighbors
#  integer, dimension(max_nibool_interfaces,num_interfaces), intent(in) :: ibool_interfaces
MESH_MPI_ARRAY_LIST = [
    ("num_interfaces", "i4"),
    ("max_nibool_interfaces", "i4"),
    ("my_neighbors", "i4"),
    ("nibool_interfaces", "i4"),
    ("ibool_interfaces", "i4"),
]

# --- data structure of proc000*_reg1_boundary_teleseismic.bin
BOUNDARY_ARRAY_LIST = [
    ("nspec2D_teleseismic_xmin", "i4"),
    ("ibelm_teleseismic_xmin", "i4"),
    ("area_teleseismic_xmin", "f4"),
    ("nspec2D_teleseismic_xmax", "i4"),
    ("ibelm_teleseismic_xmax", "i4"),
    ("area_teleseismic_xmax", "f4"),
    ("nspec2D_teleseismic_ymin", "i4"),
    ("ibelm_teleseismic_ymin", "i4"),
    ("area_teleseismic_ymin", "f4"),
    ("nspec2D_teleseismic_ymax", "i4"),
    ("ibelm_teleseismic_ymax", "i4"),
    ("area_teleseismic_ymax", "f4"),
    ("nspec2D_teleseismic_zmin", "i4"),
    ("ibelm_teleseismic_zmin", "i4"),
    ("area_teleseismic_zmin", "f4"),
]
