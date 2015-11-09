#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Postprocess for SEM3d_globe
    1. combine event kernels

"""
import os.path
import sys
import numpy as np
from fortranfile import FortranFile

#======
# constants related to SEM setup
NGLLX = NGLLY = NGLLZ = 5
NGLLCUBE = NGLLX*NGLLY*NGLLZ
R_EARTH_M = 6371000.0
R_EARTH_KM = 6371.0

def proc_filename(pdir, iproc, iregion, tag):
    """create database file name
    """
    return "%s/proc%06d_reg%1d_%s.bin" % (pdir, iproc, iregion, tag)

def read_mesh(pdir, iproc, iregion):
    """read mesh topology
    """
    f = FortranFile(proc_filename(pdir, iproc, iregion, "solver_data")
    
    nspec = f.readInts()
    nglob = f.readInts()

    glob_x = f.readReals()
    glob_y = f.readReals()
    glob_z = f.readReals()

    ngll = NGLLCUBE * nspec
    gll_xyz = np.empty((3, ngll))

    ibool = f.readReals()
    for i in range(nglob):
        iglob = ibool[i]
        gll_xyz[0, iglob] = glob_x[iglob]
        gll_xyz[1, iglob] = glob_y[iglob]
        gll_xyz[2, iglob] = glob_z[iglob]

    f.close()

    return nspec, nglob, ngll, gll_xyz

def read_gll(pdir, iproc, iregion, tag):
    """read file of GLL points  
    """
    f = FortranFile(proc_filename(pdir, iproc, iregion, tag)
    
    gll = f.readReals()

    f.close()

    return gll

def write_gll(pdir, iproc, iregion, tag, gll):
    """read file of GLL points  
    """
    f = FortranFile(proc_filename(pdir, iproc, iregion, tag), mode='w')
    
    f.writeReals(gll, prec='f')

    f.close()

def read_cmt(fname):
    """read CMTSOLUTION file
    """
    with open(fname, 'r') as f:
        cmt = [ x.split() for x in f.readlines() 
                if not(x.startswith('#')) ]
    year  = cmt[0][1]
    month  = cmt[0][2]
    day = cmt[0][3]
    hour  = cmt[0][4]
    minute = cmt[0][5]
    second = cmt[0][6]

    event_id = cmt[1][2]
    time_shift = float(cmt[2][2])
    half_duration = float(cmt[3][2])
    lat = float(cmt[4][1])
    lon = float(cmt[5][1])
    depth = float(cmt[6][1])

    # centroid time 
    isotime = '{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z'.format(
            year, month, day, hour, minute, second)
    centroid_time = UTCDateTime(isotime) + time_shift

    # moment tensor 
    # basis: (r,theta,phi) corresponds to (up,south,east)
    Mrr = float(cmt[7][1])
    Mtt = float(cmt[8][1])
    Mpp = float(cmt[9][1])
    Mrt = float(cmt[10][1])
    Mrp = float(cmt[11][1])
    Mtp = float(cmt[12][1])
    M = [[Mrr, Mrt, Mrp], [Mrt, Mtt, Mtp], [Mrp, Mtp, Mpp]]

    gcmt = {'code': event_id,
            'centroid_time':str(centroid_time),
            'half_duration':half_duration,
            'latitude':lat, 'longitude':lon, 'depth':depth, 
            'moment_tensor':M }

    return gcmt

#===== combine event kernels
def combine_kernels(nproc, mesh_dir, event_dirs, kernel_tags,
    out_dir, iregion=1, mask_source=True, radius_km=200):
    """Combine event kernels.

    Data folder structure must be:
       <event_dir>/DATABASES_MPI/proc*_reg?_solver_data.bin
                                /proc*_reg?_<kernel_tag>.bin
                  /DATA/CMTSOLUTION

    Parameters
    ----------
    nproc: int
        total number of mesh slices
    iregion: int {1, 2, 3}
        mesh region index as specified in SEM3d_globe
    sigma_km: scalar, optional
        source mask width specified in mask = 1 - exp(-dist**2/2/sigma**2)  

    """
    # get source locations if mask_source=True
    nevt = len(event_dirs)
    if mask_source:
        # initialize pyproj objects
        import pyproj
        geod = pyproj.Geod(ellps='WGS84')
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        # event x,y,z
        event_xyz = np.zeros((3, nevt))
        for idx, event_dir in enumerate(event_dirs):
            gcmt = read_cmt("%s/DATA/CMTSOLUTION" % (event_dir))
            evlat = gcmt['latitude']
            evlon = gcmt['longitude']
            evalt = -1000.0 * gcmt['depth']
            # lla to ECEF (meters)
            evx, evy, evz = pyproj.transform(lla, ecef, evlon, evlat, evalt)
            event_xyz[:, idx] = np.array([evx, evy, evz]) / R_EARTH_M
        # mask Gaussian width
        two_sigma_sq = 2.0 * (raidus_km/R_EARTH_KM)**2

    # process each slice
    mask = 1.0
    for iproc in range(nproc):
        # read mesh
        nspec, nglob, ngll, gll_xyz = read_mesh(mesh_dir, iproc, iregion)
        # process each kernel type
        for tag in kernel_tags:
            gll_sum = np.zeros(ngll)
            # read each event kernel
            for idx, event_dir in enumerate(event_dirs):
                # creat mask
                if mask_source:
                    # distances from all gll points to the source location
                    dist_sq = np.sum( 
                            (gll_xyz - event_xyz[:,(idx,)*ngll])**2, axis=0)
                    mask = 1.0 - np.exp( -dist_sq/two_sigma_sq)
                    gll_sum += read_gll(event_dir, iproc, iregion, tag) * mask
            # output file
            gll_sum /= nevt
            write_gll(out_dir, iproc, iregion, tag, gll_sum)
