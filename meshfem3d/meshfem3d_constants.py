#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
parameters from constants.h
"""

#-- radius of the Earth
R_EARTH = 6371000.0
R_EARTH_KM = 6371.0

#--- GLL

#! number of GLL points in each direction of an element (degree plus one)
NGLLX = 5
NGLLY = NGLLX
NGLLZ = NGLLX

#! mid-points inside a GLL element
MIDX = int((NGLLX-1)/2)
MIDY = int((NGLLY-1)/2)
MIDZ = int((NGLLZ-1)/2)

#! for the Gauss-Lobatto-Legendre points and weights
GAUSSALPHA = 0
GAUSSBETA = 0

#--- define flag for elements
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

#--- data structure of proc000*_reg1_solver_data.bin
MESH_ARRAY_LIST = [                       
  ('nspec','i4')                      ,
  ('nglob','i4')                      ,
  ('x','f4')                          ,
  ('y','f4')                          ,
  ('z','f4')                          ,
  ('ibool','i4')                      ,
  ('idoubling','i4')                  ,
  ('ispec_is_tiso','i4')              ,
  ('DxiDx','f4')                      ,
  ('DxiDy','f4')                      ,
  ('DxiDz','f4')                      ,
  ('DetaDx','f4')                     ,
  ('DetaDy','f4')                     ,
  ('DetaDz','f4')                     ,
  ('DgammaDx','f4')                   ,
  ('DgammaDy','f4')                   ,
  ('DgammaDz','f4')                   ,
  ]

#--- data structure of proc000*_reg1_boundary_teleseismic.bin
BOUNDARY_ARRAY_LIST = [                       
  ('nspec2D_teleseismic_xmin','i4')        ,
  ('ibelm_teleseismic_xmin', 'i4')         ,
  ('area_teleseismic_xmin', 'f4')    ,

  ('nspec2D_teleseismic_xmax','i4')        ,
  ('ibelm_teleseismic_xmax', 'i4')         ,
  ('area_teleseismic_xmax', 'f4')    ,

  ('nspec2D_teleseismic_ymin','i4')        ,
  ('ibelm_teleseismic_ymin', 'i4')         ,
  ('area_teleseismic_ymin', 'f4')    ,

  ('nspec2D_teleseismic_ymax','i4')        ,
  ('ibelm_teleseismic_ymax', 'i4')         ,
  ('area_teleseismic_ymax', 'f4')    ,

  ('nspec2D_teleseismic_zmin','i4')      ,
  ('ibelm_teleseismic_zmin', 'i4')       ,
  ('area_teleseismic_zmin', 'f4')  ,
  ]
