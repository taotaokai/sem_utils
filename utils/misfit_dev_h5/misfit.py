#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import warnings

# import os.path
import re
import datetime
from collections import OrderedDict
import multiprocessing
import queue

#
import numpy as np
import scipy
import scipy.signal as signal

# from scipy import interpolate
#
import tables as pt
import pandas as pd

# import hashlib

#
from obspy import UTCDateTime, read, Trace, Stream

# from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
from obspy.taup import TauPyModel
from obspy.imaging.beachball import beach
from obspy.signal.invsim import cosine_sac_taper
from obspy.signal.interpolation import lanczos_interpolation

#
import pyproj

# from lanczos_interp1 import lanczos_interp1

# from obspy.signal.interpolation import lanczos_interpolation

#
import matplotlib

matplotlib.use("pdf")
from matplotlib import colors, ticker, cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import cartopy.crs as ccrs
import cartopy.feature as cfeature

#
import yaml

# NOTE
# 1. spectrum relation between DFT and FT
#   x(n*dt): DFT[x]*dt ~ FT[x], IDFT[FT[x]]/dt ~ x

_DEBUG = False
# _DEBUG = True

#
# ====== utility functions
#


def is_equal(lst):
    return len(lst) >= 2 and [lst[0]] * len(lst) == lst


def stf_gauss_spectrum(f, tau):
    """
    spectrum of the Gaussian STF of unit area:
      stf(t;t0,tau) = 1/sqrt(PI)/tau * exp(-((t-t0)/tau)^2)
      here, we take t0 = 0, then
      F_stf = exp(- pi^2 * f^2 * tau^2)
    """
    F_src = np.exp(-np.pi**2 * f**2 * tau**2)
    # F_ds_dt0 = -2.0j * np.pi * f * F_src
    # F_ds_dtau = -2.0 * (np.pi*f)**2 * tau * F_src
    # return F_src, F_ds_dt0, F_ds_dtau
    return F_src


def stf_gauss_spectrum_der(f, tau):
    """
    spectrum of the derivatives of Gaussian STF
    """
    F_src = stf_gauss_spectrum(f, tau)
    F_ds_dt0 = -2.0j * np.pi * f * F_src
    F_ds_dtau = -2.0 * (np.pi * f) ** 2 * tau * F_src
    return F_ds_dt0, F_ds_dtau


def cosine_taper(x, xc):
    """cosine taper at two ends stop,pass, pass,stop
    xc: (array-like)
        stop,pass[,pass,stop]
    x: scalar or array-like
        sample points
    """
    nc = len(xc)
    if np.isscalar(x):
        x = np.array(
            [
                x,
            ]
        )
    else:
        x = np.array(x)
    y = np.ones(len(x))

    if nc == 2:  # sided taper
        l = xc[1] - xc[0]  # taper width
        if l == 0:
            raise ValueError("pass and stop values cannot be the same.")
        elif l > 0:  # tapered at left side
            idx = x < xc[0]
            y[idx] = 0.0
            idx = (xc[0] <= x) & (x <= xc[1])
            y[idx] = 0.5 - 0.5 * np.cos(np.pi * (x[idx] - xc[0]) / l)
            idx = x > xc[1]
            y[idx] = 1.0
        else:  # tapered at right side
            idx = x > xc[0]
            y[idx] = 0.0
            idx = (xc[1] <= x) & (x <= xc[0])
            y[idx] = 0.5 + 0.5 * np.cos(np.pi * (x[idx] - xc[1]) / l)
            idx = x < xc[1]
            y[idx] = 1.0
    elif nc == 4:  # two sided taper
        if not (xc[0] < xc[1] < xc[2] < xc[3]):
            raise ValueError("4 cutoff values must be in increasing order.")
        else:
            idx = x <= xc[0]
            y[idx] = 0.0

            idx = (xc[0] < x) & (x < xc[1])
            y[idx] = 0.5 - 0.5 * np.cos(np.pi * (x[idx] - xc[0]) / (xc[1] - xc[0]))

            idx = (xc[1] <= x) & (x <= xc[2])
            y[idx] = 1.0

            idx = (xc[2] < x) & (x < xc[3])
            y[idx] = 0.5 + 0.5 * np.cos(np.pi * (x[idx] - xc[2]) / (xc[3] - xc[2]))

            idx = x > xc[3]
            y[idx] = 0.0
    else:
        raise ValueError("number of cutoff values must be either 2 or 4.")

    # restore return value to scalar when input x is a scalar
    if len(y) == 1:
        y = y[0]

    return y


def spheredist(lat0, lon0, lat1, lon1):
    d2r = np.pi / 180.0
    lat0, lon0 = lat0 * d2r, lon0 * d2r
    lat1, lon1 = lat1 * d2r, lon1 * d2r
    # v0 = [np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
    # v1 = [np.cos(lat1)*np.cos(lon1), np.cos(lat1)*np.sin(lon1), np.sin(lat1)]
    return np.arccos(
        np.sin(lat0) * np.sin(lat1) + np.cos(lat0) * np.cos(lat1) * np.cos(lon1 - lon0)
    )


def centerMap(lats, lons, scale):
    """
    https://stackoverflow.com/questions/13240642/automatically-center-matplotlib-basemap-onto-data
    """
    # NOTE: only works for northern hemisphere!
    # Assumes -90 < Lat < 90 and -180 < Lon < 180, and
    # latitude and logitude are in decimal degrees
    earthRadius = 6378100.0  # earth's radius in meters
    northLat = max(lats)
    southLat = min(lats)
    westLon = max(lons)
    eastLon = min(lons)
    # average between max and min longitude
    lon0 = ((westLon - eastLon) / 2.0) + eastLon
    # lon0 = np.median(lons)
    # a = the height of the map
    b = spheredist(northLat, westLon, northLat, eastLon) * earthRadius / 2
    c = spheredist(northLat, westLon, southLat, lon0) * earthRadius
    # use pythagorean theorom to determine height of plot
    mapH = pow(pow(c, 2) - pow(b, 2), 1.0 / 2)
    arcCenter = (mapH / 2) / earthRadius
    lat0 = southLat + arcCenter * 180.0 / np.pi
    # distance between max E and W longitude at most souther latitude
    # widest part on the map, either at the equator or at the latitude closest to the equator
    minLat = min(abs(southLat), abs(northLat))
    if np.sign(southLat) != np.sign(northLat):
        minLat = 0
    mapW = spheredist(minLat, westLon, minLat, eastLon) * earthRadius
    return lat0, lon0, mapW * scale, mapH * scale


# def vector2latlon_deg(v):
#     lon = np.arctan2(v[1], v[0])
#     latc = np.arctan2(v[2], (v[0] ** 2 + v[1] ** 2) ** 0.5)
#     f = 1 / 299.8
#     factor = 1.0 / (1 - f) ** 2
#     lat = np.arctan(factor * np.tan(latc))
#     return np.rad2deg(lat), np.rad2deg(lon)


# def geodetic_lat2geocentric_lat(geodetic_lat):
#     f = 1 / 299.8
#     factor = (1 - f) ** 2
#     return np.arctan(factor * np.tan(geodetic_lat))
#
#
# def get_dist_from_mesh_boundary(lat_center, lon_center, gamma_rot, lat_test, lon_test):
#     """get mesh specific xi/eta for (lat_test,lon_test)
#     useful for testing whether a point lies inside the SEM mesh region
#     """
#     lat0 = np.deg2rad(lat_center)
#     lon0 = np.deg2rad(lon_center)
#     # assume zero altitude from the reference ellipsoid
#     theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat0)
#     phi = lon0
#     # radial/easting/northing direction at (lat_center, lon_center)
#     v0_r = np.array(
#         [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
#     )
#     v0_e = np.array([-np.sin(phi), np.cos(phi), 0])
#     v0_n = np.array(
#         [-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)]
#     )
#
#     # rotate (v0_e, v0_n) to (v0_xi, v0_eta) through v0_r by gamma_rot counter-clockwise
#     gamma = np.deg2rad(gamma_rot)
#     v0_xi = np.cos(gamma) * v0_e + np.sin(gamma) * v0_n
#     v0_eta = -np.sin(gamma) * v0_e + np.cos(gamma) * v0_n
#     # print(v_xi, v_eta)
#
#     # test point
#     lat1 = np.deg2rad(lat_test)
#     lon1 = np.deg2rad(lon_test)
#     theta = 0.5 * np.pi - geodetic_lat2geocentric_lat(lat1)
#     phi = lon1
#     v1 = np.array(
#         [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
#     )
#     # project to v0_r, v0_xi, v0_eta
#     l_r = np.dot(v0_r, v1)
#     l_xi = np.dot(v0_xi, v1)
#     l_eta = np.dot(v0_eta, v1)
#     angle_xi = np.rad2deg(np.arctan2(l_xi, l_r))
#     angle_eta = np.rad2deg(np.arctan2(l_eta, l_r))
#
#     return angle_xi, angle_eta


#
# ======== for use internally
#


class Channel(pt.IsDescription):
    network = pt.StringCol(10, pos=0)
    station = pt.StringCol(10, pos=1)
    location = pt.StringCol(2, pos=2)
    channel = pt.StringCol(3, pos=3)
    latitude = pt.Float64Col(pos=4)
    longitude = pt.Float64Col(pos=5)
    elevation = pt.Float64Col(pos=6)
    depth = pt.Float64Col(pos=7)
    azimuth = pt.Float64Col(pos=8)
    dip = pt.Float64Col(pos=9)
    starttime = pt.StringCol(30, pos=10)
    endtime = pt.StringCol(30, pos=11)


class Window(pt.IsDescription):
    network = pt.StringCol(10, pos=0)
    station = pt.StringCol(10, pos=1)
    type = pt.StringCol(10, pos=2)
    phase = pt.StringCol(20, pos=3)
    cmpnm = pt.StringCol(1, pos=4)
    starttime = pt.StringCol(30, pos=5)
    endtime = pt.StringCol(30, pos=6)
    taper = pt.Float64Col(pos=7)
    butter_N = pt.Int32Col(pos=8)
    butter_Wn = pt.Float64Col(2, pos=9)
    cmpaz = pt.Float64Col(pos=10)
    cmpdip = pt.Float64Col(pos=11)
    weight = pt.Float64Col(pos=12)
    cc0 = pt.Float64Col(pos=13)
    cc_time_shift = pt.Float64Col(
        pos=14
    )  # syn(t-dt): positive dt means shiftting syn right
    cc_max = pt.Float64Col(pos=15)  # cc after time shift
    amp_ratio_cc0 = pt.Float64Col(
        pos=16
    )  # (w*obs, w*syn)/(w*obs, w*obs), (*,*): inner product
    amp_ratio_ccmax = pt.Float64Col(pos=17)  # after time shift
    SNR = pt.Float64Col(pos=18)
    noise_maxamp = pt.Float64Col(pos=19)
    obs_maxamp = pt.Float64Col(pos=20)
    syn_maxamp = pt.Float64Col(pos=21)
    id = pt.StringCol(128, pos=22)
    valid = pt.BoolCol(pos=23)
    proj_matrix = pt.Float64Col(
        shape=(3, 3), pos=24
    )  # from ENZ to the measured component


class Source(pt.IsDescription):
    header = pt.StringCol(
        100, pos=0
    )  # e.g. PDEW 2010 06 18 23 09 34.3000 13.20 93.09 20.0 6.1 5.9 ANDAMAN ISLANDS, INDIA R
    id = pt.StringCol(30, pos=1)
    t0 = pt.StringCol(30, pos=2)  # t0
    tau = pt.Float64Col(pos=3)  # tau
    latitude = pt.Float64Col(pos=4)
    longitude = pt.Float64Col(pos=5)
    depth = pt.Float64Col(pos=6)
    mt_rtp = pt.Float64Col(shape=(3, 3), pos=7)  # in CMT r,theta,phi coordinates
    xs = pt.Float64Col(shape=(3), pos=8)  # xs in ECEF coordinate
    mt = pt.Float64Col(shape=(3, 3), pos=9)  # in ECEF coordinate
    # force_vector =
    # gradient
    grad_xs = pt.Float64Col(shape=(3), pos=11)  # dchi_dxs
    grad_mt = pt.Float64Col(shape=(3, 3), pos=12)  # dchi_dmt
    # search direction for xs, mt
    dxs = pt.Float64Col(shape=(3), pos=13)
    dmt = pt.Float64Col(shape=(3, 3), pos=14)
    # dt0 = dtau = 1
    # optimal step length
    step_xs = pt.Float64Col(pos=15)
    step_mt = pt.Float64Col(pos=16)
    step_t0 = pt.Float64Col(pos=17)
    step_tau = pt.Float64Col(pos=18)
    # NOTE the update source parameters:
    #   xs <- xs + step_xs * dxs
    #   mt <- mt + step_mt * dmt
    #   t0 <- t0 + step_t0
    #   tau <- tau + step_tau


def _get_obs_ENZ(g_sta, obs_tag):
    """get observed seismogram in array(ENZ,nt)"""
    if obs_tag not in g_sta:
        msg = f"{g_sta._v_name}: {obs_tag} does not exist, skip"
        raise KeyError(msg)
    obs = g_sta[obs_tag]

    # check data type
    obs_type = obs.attrs["type"]
    if obs_type != "VEL" and obs_type != "DISP":
        raise AssertionError(f"wrong types: obs({obs_type}), either VEL or DISP")

    # check channels' info
    obs_channels = obs.attrs["channels"]
    obs_Zchan = [cha for cha in obs_channels if cha["name"].decode()[-1] == "Z"]
    obs_Hchan = [cha for cha in obs_channels if cha["name"].decode()[-1] != "Z"]
    if len(obs_Zchan) != 1 or obs_Zchan[0]["dip"] != -90:
        msg = f"{g_sta._v_name} has problematic Z channel info: {obs_Zchan}, skip"
        raise AssertionError(msg)
    no_Hchan = False
    if len(obs_Hchan) == 0:
        msg = f"{g_sta._v_name} has no horizontal channels"
        warnings.warn(msg)
        no_Hchan = True
    else:
        assert len(obs_Hchan) == 2

    # rotate obs to ENZ
    # projection matrix: obs = proj * ENZ => ENZ = inv(proj) * obs
    data_nt = obs.attrs["npts"]
    obs_ENZ = np.zeros((3, data_nt), dtype=float)
    if no_Hchan:
        obs_ENZ[0:2, :] = 0.0
        obs_ENZ[2, :] = obs[0, :]  # only Z component
    else:
        assert obs.shape == obs_ENZ.shape
        proj_matrix = np.zeros((3, 3))
        for i in range(3):
            chan = obs_channels[i]
            sin_az = np.sin(np.deg2rad(chan["azimuth"]))
            cos_az = np.cos(np.deg2rad(chan["azimuth"]))
            sin_dip = np.sin(np.deg2rad(chan["dip"]))
            cos_dip = np.cos(np.deg2rad(chan["dip"]))
            # column vector = obs channel polarization
            proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
            proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
            proj_matrix[i, 2] = -sin_dip  # proj to Z (Up)
        # inverse projection matrix: ENZ = inv(proj) * obs
        inv_proj = np.linalg.inv(proj_matrix)
        obs_ENZ = np.dot(inv_proj, obs)

    return obs_ENZ, no_Hchan


def _get_syn_ENZ(g_sta, syn_tag):

    if syn_tag not in g_sta:
        msg = f"{g_sta._v_name}: {syn_tag} does not exist, skip"
        raise KeyError(msg)
    syn = g_sta[syn_tag]

    # check data type
    syn_type = syn.attrs["type"]
    if syn_type != "DISP":
        raise AssertionError(f"wrong types: syn({syn_type}) must be DISP")

    # check channels' info
    syn_channels = syn.attrs["channels"]
    assert len(syn_channels) == 3

    # rotate syn to ENZ
    data_nt = syn.attrs["npts"]
    data_fs = syn.attrs["sampling_rate"]
    data_dt = 1.0 / data_fs
    syn_ENZ = np.zeros((3, data_nt))
    syn_ENZ[:] = syn[:]
    proj_matrix = np.zeros((3, 3))
    for i in range(3):
        chan = syn_channels[i]
        sin_az = np.sin(np.deg2rad(chan["azimuth"]))
        cos_az = np.cos(np.deg2rad(chan["azimuth"]))
        sin_dip = np.sin(np.deg2rad(chan["dip"]))
        cos_dip = np.cos(np.deg2rad(chan["dip"]))
        # column vector = obs channel polarization
        proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
        proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
        proj_matrix[i, 2] = -sin_dip  # proj to Z
    # inverse projection matrix: ENZ = inv(proj) * obs
    inv_proj = np.linalg.inv(proj_matrix)
    syn_ENZ[:] = np.dot(inv_proj, syn_ENZ)

    # # apply time derivative (when data is vel)
    # if to_vel:
    #     npad = 100  # t_pad = 100 / fs, response drop to ~1% peak value
    #     nfft = scipy.fft.next_fast_len(data_nt + npad)
    #     freqs = np.fft.rfftfreq(nfft, d=data_dt)
    #     Fsyn *= 2j * np.pi * freqs
    #     syn_ENZ[:] = np.fft.irfft(Fsyn, nfft)[:, :data_nt]

    return syn_ENZ


def _extract_obs_syn_ENZ(
    obs_data,
    obs_attrs,
    syn_data,
    syn_attrs,  # , event_tau
    # sta_attrs, obs_data, obs_attrs, syn_data, syn_attrs  # , event_tau
):
    # net = sta_attrs["network"]
    # sta = sta_attrs["station"]
    # stnm = f"{net}_{sta}"

    # check data type
    obs_type = obs_attrs["type"]
    syn_type = syn_attrs["type"]
    case1 = obs_type == "VEL" and syn_type == "DISP"
    case2 = obs_type == "DISP" and syn_type == "DISP"
    if not (case1 or case2):
        raise AssertionError(f"wrong types: obs({obs_type}), syn({syn_type})")

    # check channels' info
    obs_channels = obs_attrs["channels"]
    syn_channels = syn_attrs["channels"]
    obs_Zchan = [cha for cha in obs_channels if cha["name"].decode()[-1] == "Z"]
    obs_Hchan = [cha for cha in obs_channels if cha["name"].decode()[-1] != "Z"]
    if len(obs_Zchan) != 1 or obs_Zchan[0]["dip"] != -90:
        # msg = f"{stnm} has problematic Z channel info: {obs_Zchan}, skip"
        msg = f"problematic Z channel info: {obs_Zchan}, skip"
        raise AssertionError(msg)
    no_Hchan = False
    if len(obs_Hchan) == 0:
        # msg = f"{stnm} has no horizontal channels"
        msg = f"no horizontal channels"
        warnings.warn(msg)
        no_Hchan = True
    else:
        assert len(obs_Hchan) == 2
    assert len(syn_channels) == 3

    # data_tb = obs.attrs["starttime"]
    data_fs = obs_attrs["sampling_rate"]
    data_dt = 1.0 / data_fs
    data_nt = obs_attrs["npts"]
    # data_te = data_tb + (data_nt - 1) * data_dt
    # data_times = np.arange(data_nt, dtype=float) * data_dt
    # noise_te = (first_arrtime - data_tb) - cfg_noise_before_first_arrival
    # noise_idx1 = int(noise_te * data_fs)  # 0:idx1 as noise

    # rotate obs to ENZ
    # projection matrix: obs = proj * ENZ => ENZ = inv(proj) * obs
    obs_ENZ = np.zeros((3, data_nt))
    if no_Hchan:
        obs_ENZ[2, :] = obs_data[0, :]  # only Z component
    else:
        assert obs_data.shape == obs_ENZ.shape
        proj_matrix = np.zeros((3, 3))
        for i in range(3):
            chan = obs_channels[i]
            sin_az = np.sin(np.deg2rad(chan["azimuth"]))
            cos_az = np.cos(np.deg2rad(chan["azimuth"]))
            sin_dip = np.sin(np.deg2rad(chan["dip"]))
            cos_dip = np.cos(np.deg2rad(chan["dip"]))
            # column vector = obs channel polarization
            proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
            proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
            proj_matrix[i, 2] = -sin_dip  # proj to Z (Up)
        # inverse projection matrix: ENZ = inv(proj) * obs
        inv_proj = np.linalg.inv(proj_matrix)
        obs_ENZ = np.dot(inv_proj, obs_data)

    # rotate syn to ENZ
    syn_ENZ = np.zeros((3, data_nt))
    syn_ENZ[:] = syn_data[:]
    proj_matrix = np.zeros((3, 3))
    for i in range(3):
        chan = syn_channels[i]
        sin_az = np.sin(np.deg2rad(chan["azimuth"]))
        cos_az = np.cos(np.deg2rad(chan["azimuth"]))
        sin_dip = np.sin(np.deg2rad(chan["dip"]))
        cos_dip = np.cos(np.deg2rad(chan["dip"]))
        # column vector = obs channel polarization
        proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
        proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
        proj_matrix[i, 2] = -sin_dip  # proj to Z
    # inverse projection matrix: ENZ = inv(proj) * obs
    inv_proj = np.linalg.inv(proj_matrix)
    syn_ENZ[:] = np.dot(inv_proj, syn_ENZ)

    # np.savez('extract', syn_ENZ=syn_ENZ)

    # # apply source time function and/or time derivative
    # if syn_attrs["is_grn"] or obs_type == "VEL":
    #     npad = int(5 * event_tau * data_fs)
    #     # npad = data_nt  #DEBUG
    #     nfft = scipy.fft.next_fast_len(data_nt + npad)
    #     freqs = np.fft.rfftfreq(nfft, d=data_dt)
    #     Fsyn = np.fft.rfft(syn_ENZ, nfft)
    #     if syn_attrs["is_grn"]:
    #         # source spectrum (moment-rate function)
    #         # print("event_tau=", event_tau)
    #         F_src = stf_gauss_spectrum(freqs, event_tau)
    #         Fsyn *= F_src
    #     if obs_type == "VEL":
    #         Fsyn *= 2j * np.pi * freqs
    #     syn_ENZ[:] = np.fft.irfft(Fsyn, nfft)[:, :data_nt]

    return obs_ENZ, syn_ENZ, no_Hchan


def _measure_adj_one_sta(
    sta_attrs,
    obs_data,
    obs_attrs,
    syn_data,
    syn_attrs,
    wins,
    event_t0,
    event_tau,
    cfg_noise_before_first_arrival,
    cci_dt,
    misfit_type,
    weight_param,
    adj_band_code,
):

    net = sta_attrs["network"]
    sta = sta_attrs["station"]
    stnm = f"{net}_{sta}"
    # loc = attrs['location']
    first_arrtime = sta_attrs["first_arrtime"]

    try:
        obs_ENZ, syn_ENZ, no_Hchan = _extract_obs_syn_ENZ(
            obs_data,
            obs_attrs,
            syn_data,
            syn_attrs,  # , event_tau
            # sta_attrs, obs_data, obs_attrs, syn_data, syn_attrs  # , event_tau
        )
        # debug
        # np.savez(f'{net}_{sta}_new.npz', obs_ENZ=obs_ENZ, syn_ENZ=syn_ENZ, no_Hchan=no_Hchan)
    except Exception as e:
        msg = f"{stnm}: failed to get obs_ENZ,syn_ENZ ({e})"
        warnings.warn(msg)
        return None

    data_tb = obs_attrs["starttime"]
    data_fs = obs_attrs["sampling_rate"]
    data_dt = 1.0 / data_fs
    data_nt = obs_attrs["npts"]
    data_te = data_tb + (data_nt - 1) * data_dt
    data_times = np.arange(data_nt, dtype=float) * data_dt
    noise_te = (first_arrtime - data_tb) - cfg_noise_before_first_arrival
    # noise_idx1 = int(noise_te * data_fs)  # 0:idx1 as noise

    # pre-filter parameters
    pre_filters = obs_attrs["filter"]
    # data_butter_N = obs.attrs["filter"]["N"]
    # data_butter_Wn = obs.attrs["filter"]["Wn"]

    obs_type = obs_attrs["type"]
    syn_type = syn_attrs["type"]
    if obs_type == "VEL" and syn_type == "DISP":
        syn_to_vel = True
    elif obs_type == "DISP" and syn_type == "DISP":
        syn_to_vel = False
    else:
        msg = f"{stnm}: unknown obs/syn types ({obs_type}/{syn_type})"
        raise ValueError(msg)

    conv_stf = False
    if syn_attrs["is_grn"]:
        conv_stf = True
        assert syn_attrs["origin_time"] == event_t0

    # sum of adjoint sources from all windows
    dchi_du = np.zeros((3, data_nt))

    for win in wins:
        # component
        cmpnm = win["cmpnm"].decode()
        if no_Hchan and cmpnm != "Z":
            warnings.warn("skip non-vertical window since only vertical data present")
            win["valid"] = False
            win.update()
            continue
        cmpaz = win["cmpaz"]
        cmpdip = win["cmpdip"]

        # filter
        butter_N = win["butter_N"]
        butter_Wn = win["butter_Wn"]
        # fft
        pad_length = max(20.0, 2.0 / min(butter_Wn)) # pad both ends before apply filter to reduce edge effect
        npad = int(pad_length * data_fs)
        # nfft = scipy.fft.next_fast_len(data_nt + npad)
        nfft = scipy.fft.next_fast_len(data_nt + 2 * npad)
        taper = cosine_taper(np.arange(nfft), [0, npad, npad + data_nt, nfft-1]) # applied before filtering
        freqs = np.fft.rfftfreq(nfft, d=data_dt)
        # window filter
        sos = scipy.signal.butter(
            butter_N, butter_Wn, "bandpass", fs=data_fs, output="sos"
        )
        _, filter_h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
        Fw = abs(filter_h)  # zero-phase filter response
        # pre-filter applied to obs & syn
        filter_h = np.ones_like(freqs)
        for filter in pre_filters:
            sos = scipy.signal.butter(
                filter["N"],
                filter["Wn"],
                filter["type"],
                fs=data_fs,
                output="sos",
            )
            _, h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
            filter_h *= abs(h)
        Fd = abs(filter_h)
        # time derivative in frequency domain
        Ft = 2j * np.pi * freqs
        # time window
        win_tb = UTCDateTime(win["starttime"])
        win_te = UTCDateTime(win["endtime"])
        win_taper = win["taper"]

        # window function
        b = win_tb - data_tb
        e = win_te - data_tb
        win_len = e - b
        width = win_len * min(win_taper, 0.4)
        win_c = [b, b + width, e - width, e]
        win_func = cosine_sac_taper(data_times, win_c)
        # win_idx = (data_times >= b + width) & (data_times <= e - width)

        # component
        if cmpnm in ["Z", "R", "T"]:
            sin_az = np.sin(np.deg2rad(cmpaz))
            cos_az = np.cos(np.deg2rad(cmpaz))
            sin_dip = np.sin(np.deg2rad(cmpdip))
            cos_dip = np.cos(np.deg2rad(cmpdip))
            n = np.array(
                [
                    [cos_dip * sin_az],  # cos(E, comp)
                    [cos_dip * cos_az],  # cos(N, comp)
                    [-sin_dip],  # cos(Z, comp)
                ]
            )
            proj_matrix = np.dot(n, n.transpose())
        elif cmpnm == "H":  # horizontal vector 2d
            proj_matrix = np.identity(3)
            proj_matrix[2, 2] = 0.0  # zero Z component
        elif cmpnm == "F":  # full 3d vector
            proj_matrix = np.identity(3)
        else:
            msg = f"unrecognized component code ({cmpnm}), SKIP"
            warnings.warn(msg)
            win["valid"] = False
            continue

        # Fd: pre-filter used during preprocessing (e.g. rmresp)
        # Fw: bandpass filter used on this window
        # w: window function
        # S: source spectrum
        # Ft: time derivative in frequency domain
        #
        # obs = w*(Fw*Fd*d), syn = w*(Fw*Fd*u)

        # apply filter to obs, syn
        # obs = Fw * (Fd * d), obs_ENZ = Fd * d, d = disp. or vel.
        # obs_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ, nfft), nfft)[:, :data_nt]
        obs_ENZ_pad = np.pad(obs_ENZ, ((0,0), (npad, nfft-data_nt-npad)), "edge") * taper[None, :]
        obs_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ_pad, nfft), nfft)[:, npad:(npad+data_nt)]
        # syn = Fw * (Fd * [Ft] * [S] * u), (syn_ENZ = Fd * u)
        # f_syn = np.fft.rfft(syn_ENZ, nfft)
        syn_ENZ_pad = np.pad(syn_ENZ, ((0,0), (npad, nfft-data_nt-npad)), "edge") * taper[None, :]
        f_syn = np.fft.rfft(syn_ENZ_pad, nfft)
        if syn_to_vel:
            f_syn *= Ft
        if conv_stf:
            f_src = stf_gauss_spectrum(freqs, event_tau)
            f_syn *= f_src
        # syn_filt = np.fft.irfft(Fw * f_syn, nfft)[:, :data_nt]
        syn_filt = np.fft.irfft(Fw * f_syn, nfft)[:, npad:(npad+data_nt)]
        # noise
        taper_width = 0.5 / min(butter_Wn)
        noise_win = cosine_sac_taper(
            data_times, [0, taper_width, noise_te - taper_width, noise_te]
        )
        noise_idx = (data_times >= taper_width) & (data_times <= noise_te - taper_width)
        # # cut noise window from obs_ENZ before applying bandpass filter Fw
        # noise_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ * noise_win, nfft), nfft)[
        #     :, :data_nt
        # ]

        # apply window and projection
        # obs: w*(Fw*Fd*d)
        obs_filt_win = np.dot(proj_matrix, obs_filt) * win_func
        # syn: w*(Fw*Fd*u)
        syn_filt_win = np.dot(proj_matrix, syn_filt) * win_func
        # noise
        # noise_filt_win = np.dot(proj_matrix, noise_filt)[:, :noise_idx1]
        noise_filt_win = np.dot(proj_matrix, obs_filt) * noise_win

        # ------ measure SNR (based on maximum amplitude)
        Amax_obs = np.sqrt(np.max(np.sum(obs_filt_win**2, axis=0)))
        Amax_syn = np.sqrt(np.max(np.sum(syn_filt_win**2, axis=0)))
        # Amax_noise = np.sqrt(np.max(np.sum(noise_filt_win**2, axis=0)))
        Amax_noise = np.sqrt(np.max(np.sum(noise_filt_win[:, noise_idx] ** 2, axis=0)))

        if Amax_obs == 0:  # bad record
            msg = f"empty obs trace ({win}), skip"
            warnings.warn(msg)
            win["valid"] = False
            continue
        if Amax_noise == 0:
            # could occure when using synthetic seismograms as data for pointspread test
            msg = f"noise amplitude is zero for ({win}). SNR is set to 9999.0"
            warnings.warn(msg)
            snr = 9999.0
            # win["valid"] = False
            # continue
        else:
            snr = 20.0 * np.log10(Amax_obs / Amax_noise)

        # ------ measure CC time shift (between w*F*d and w*F*u)
        obs_norm = np.sqrt(np.sum(obs_filt_win**2))
        syn_norm = np.sqrt(np.sum(syn_filt_win**2))
        # window normalization factor (without dt)
        Nw = obs_norm * syn_norm
        # NOTE the order (obs,syn) is important. The positive time on
        # CC means shifting syn in the positive time direction to match
        # the observed obs, and vice verser.
        # [-(nt-1), nt) * dt
        cc = np.zeros(2 * data_nt - 1)
        for i in range(3):
            cc += signal.fftconvolve(obs_filt_win[i, :], syn_filt_win[i, ::-1], "full")
        cc /= Nw
        # -- zero-lag cc coeff.
        CC0 = cc[data_nt - 1]  # the n-th point corresponds to zero lag time
        AR0 = CC0 * syn_norm / obs_norm  # amplitude ratio syn/obs
        # DEBUG
        # print(CC0, np.sum(obs_filt_win * syn_filt_win)/obs_norm/syn_norm)
        # -- interpolate cc to finer time samples
        CC_shift_range = win_len / 2.0  # TODO: more reasonable choice?
        if data_dt < cci_dt:
            msg = f"data_dt({data_dt}) < cci_dt({cci_dt})"
            warnings.warn(msg)
            cci_dt = data_dt
        ncci = int(CC_shift_range / cci_dt)
        cci_t0 = -ncci * cci_dt
        cc_t0 = -(data_nt - 1) * data_dt  # begin time in cc
        cci = lanczos_interpolation(
            cc, cc_t0, data_dt, cci_t0, cci_dt, 2 * ncci + 1, a=20
        )
        # time shift at the maximum correlation
        imax = np.argmax(cci)
        CC_time_shift = cci_t0 + imax * cci_dt
        CCmax = cci[imax]
        ARmax = CCmax * syn_norm / obs_norm  # amplitude ratio: syn/obs

        # ------ window weighting based on SNR and misfit
        weight = win["weight"]
        if "SNR" in weight_param:
            weight *= cosine_taper(snr, weight_param["SNR"])
        if "CCmax" in weight_param:
            weight *= cosine_taper(CCmax, weight_param["CCmax"])
        if "CC0" in weight_param:
            weight *= cosine_taper(CC0, weight_param["CC0"])
        if "cc_tshift" in weight_param:
            weight *= cosine_taper(CC_time_shift, weight_param["cc_tshift"])

        # ------ measure adjoint source
        Aw = CC0 * obs_norm / syn_norm  # window amplitude raito
        if misfit_type == "cc0":
            #
            # Goodness-of-fit: (zero-lag normalized cross-correlation coefficient)
            #   cc(syn, obs) = sum(syn * obs) / sqrt(sum(syn**2) * sum(obs**2))
            #
            # obs = w * (Fw * Fd * P * d),        d = displacement or velocity.
            # syn = w * (Fw * Fd * P * u),        if d = disp.
            #       w * (Fw * Fd * Ft * P * u),   if d = vel.
            #
            #  u[3,nt]: synthetic seismogram of displacement
            #     , replaced by S * u, when u is Green's function (S: source wavelet)
            #  d[3,nt]: recorded ground displacement or velocity
            #  P[3,3]: the projection matrix (real valued), e.g. project to radial/tangential component
            #  Fd[nfft]: bandpass filter applied during pre-processing (e.g. removing instrument response, down-sampling)
            #  Fw[nfft]: bandpass filter applied for the measurment time window
            #  Ft[nfft]: time derivate (2j * pi * freqs in frequency domian)
            #  note: F * x = irfft(F * rfft(x, nfft), nfft)[:, :nt]
            #  w[nt]: time window function
            #
            # adjoint source:
            #   Dcc/Du[3,nt] = conj(Fw * Fd * [Ft]) * transpose(P) * w * (obs - Aw * syn) / obs_norm / syn_norm,
            #
            #   where,
            #           [Ft] is only applied when d is vel.
            #           obs_norm = sqrt(sum(obs**2))
            #           syn_norm = sqrt(sum(syn**2))
            #           Aw = sum(syn * obs) / sum(syn**2)
            #
            #   cc(u+du) - cc(u) = sum(Dcc/Du * du) + O(du**2)
            #
            dchiw_du = np.dot(
                np.transpose(proj_matrix),
                win_func * (obs_filt_win - Aw * syn_filt_win) / Nw,
            )
            conj_F = np.conjugate(Fw * Fd, dtype="complex")
            if obs_attrs["type"] == "VEL":
                conj_F *= np.conjugate(Ft)
            dchiw_du = np.fft.irfft(conj_F * np.fft.rfft(dchiw_du, nfft), nfft)[
                :, :data_nt
            ]
        else:
            msg = f"{stnm}: unknown misfit type ({misfit_type})"
            raise ValueError(msg)

        # add into total dchi_du
        dchi_du += weight * dchiw_du

        win["cc0"] = CC0
        win["cc_time_shift"] = CC_time_shift
        win["cc_max"] = CCmax
        win["amp_ratio_cc0"] = AR0
        win["amp_ratio_ccmax"] = ARmax
        win["noise_maxamp"] = Amax_noise
        win["obs_maxamp"] = Amax_obs
        win["syn_maxamp"] = Amax_syn
        win["SNR"] = snr
        win["weight"] = weight
        win["valid"] = True
        win["proj_matrix"] = proj_matrix

    # np.savez(f'{net}_{sta}_wins_mp.npz', wins=wins)

    # store adjoint source for this station, e.g. /NET_STA/ADJ_DISP[0:nchan, 0:npts]
    # if adj_tag in g_sta:
    #     msg = f"{adj_tag} exists in {g_sta._v_name}, overwrite!"
    #     warnings.warn(msg)
    #     h5f.remove_node(g_sta, adj_tag)
    # shape = (3, data_nt)
    # ca = h5f.create_carray(
    #     g_sta, adj_tag, pt_atom, shape, filters=pt_filters
    # )
    # ca[:] = dchi_du
    dtype = [("name", "S3"), ("azimuth", float), ("dip", float)]
    channels = [
        (f"{adj_band_code}E", 90, 0),
        (f"{adj_band_code}N", 0, 0),
        (f"{adj_band_code}Z", 0, -90),
    ]
    adj_attrs = {}
    adj_attrs["channels"] = np.array(channels, dtype=dtype)
    adj_attrs["starttime"] = data_tb
    adj_attrs["sampling_rate"] = data_fs
    adj_attrs["npts"] = data_nt
    adj_attrs["solver_tb"] = syn_attrs["solver_tb"]
    adj_attrs["solver_dt"] = syn_attrs["solver_dt"]
    adj_attrs["solver_nt"] = syn_attrs["solver_nt"]
    adj_attrs["type"] = syn_attrs["type"]
    adj_attrs["misfit"] = misfit_type

    return dchi_du, adj_attrs, wins


def _linearized_search_cc_by_stations(args):
    # h5_path, # path to misfit.h5 file
    # dm,  # dict of model values, e.g. {"dt0": None, "dtau": None, "dxs": None, "dmt": None},
    # stations,  # list of stations, e.g. ['NET.STA', 'NN,SS']
    h5_path, dm, stations = args

    # check model vectors in dm have the same length
    if (not dm) or len(dm) == 0:
        msg = "dm is empty!"
        raise ValueError(msg)
    # model_num = len(dm)

    # for model_name in dm:
    #     if model_name not in ["dt0", "dtau", "dxs", "dmt"]:
    #         msg = f"model_name[{model_name}] must be one of [t0, tau, xs, mt]"
    #         raise ValueError(msg)

    nstep = []
    for model_name in dm:
        if dm[model_name].ndim != 1:
            msg = f"dm[{model_name}] must be 1-D vector"
            raise ValueError(msg)
        nstep.append(np.size(dm[model_name]))

    if len(set(nstep)) != 1 or nstep[0] < 1:
        msg = "vectors in dm must have the same length"
        raise Exception(msg)
    nstep = nstep[0]

    h5f = pt.open_file(h5_path, "r")

    # check parameters
    if "/source" not in h5f:
        msg = '"/source" not existing, run read_cmtsolution first!'
        raise KeyError(msg)
    tbl_src = h5f.get_node("/source")
    if tbl_src.nrows == 0:
        msg = "no source information"
        raise Exception(msg)
    event = tbl_src[0]
    event_t0 = UTCDateTime(event["t0"])

    if "dtau" in dm:
        tau = dm["dtau"] + event["tau"]
        if any(tau < 0):
            msg = f"too small dm['dtau'] ({dm['dtau']}) for event['tau']={event['tau']}"
            raise ValueError(msg)

    # syn_components = ["E", "N", "Z"]

    config = h5f.root._v_attrs["config"]
    obs_tag = config["data"]["tag"]
    syn_tag = config["syn"]["tag"]
    dsyn_grp = config["syn"]["dsyn_grp"]

    # sum of weighted normalized zero-lag cc at each model grid
    wcc_sum = np.zeros(nstep)
    # sum of all windows' weighting
    weight_sum = 0.0

    if "/waveform" not in h5f:
        msg = '"/waveform" not existing, run read_data_h5 first!'
        raise KeyError(msg)
    g_wav = h5f.get_node("/waveform")

    if "/window" not in h5f:
        msg = '"/window" not existing, run setup_window/measure_adj first!'
        raise KeyError(msg)
    tbl_win = h5f.get_node("/window")

    dsyn_tags = {}
    for model_name in dm:
        if model_name not in ["dt0", "dtau"]:
            dsyn_tags[model_name] = f"{dsyn_grp}/{model_name}"

    # print(dsyn_tags)

    for g_sta in g_wav:
        # print(f"cc_linearize for {g_sta._v_name}")
        if g_sta._v_name not in stations:
            continue

        attrs = g_sta._v_attrs
        net = attrs["network"]
        sta = attrs["station"]

        windows = tbl_win.read_where(
            f'(network == b"{net}") & (station == b"{sta}") & (valid == True) & (weight > 0)'
        )
        if len(windows) == 0:
            msg = f"No windows found for {g_sta._v_name}, SKIP"
            warnings.warn(msg)
            continue

        assert obs_tag in g_sta
        assert syn_tag in g_sta
        for _, tag in dsyn_tags.items():
            assert tag in g_sta

        obs = g_sta[obs_tag]
        syn = g_sta[syn_tag]

        obs_type = obs.attrs["type"]
        syn_type = syn.attrs["type"]
        if obs_type == "VEL" and syn_type == "DISP":
            syn_to_vel = True
        elif obs_type == "DISP" and syn_type == "DISP":
            syn_to_vel = False
        else:
            msg = f"unknown obs/syn types ({obs_type}/{syn_type}) in {g_sta._v_name}"
            raise ValueError(msg)

        for _, dsyn_tag in dsyn_tags.items():
            assert g_sta[dsyn_tag].attrs["is_grn"] == g_sta[syn_tag].attrs["is_grn"]
            assert g_sta[dsyn_tag].attrs["type"] == g_sta[syn_tag].attrs["type"]

        try:
            dsyn_dm = {}
            obs_ENZ, no_Hchan = _get_obs_ENZ(g_sta, obs_tag)
            syn_ENZ = _get_syn_ENZ(g_sta, syn_tag)
            for m, tag in dsyn_tags.items():
                dsyn_dm[m] = _get_syn_ENZ(g_sta, tag)
        except Exception as e:
            msg = f"failed to get {obs_tag},{syn_tag},{dsyn_tags} from {g_sta._v_name}, ({e})"
            warnings.warn(msg)
            continue

        if "dtau" in dm:
            assert g_sta[syn_tag].attrs["is_grn"] == True
            # for m, tag in dsyn_tags:
            #     assert g_sta[tag].attrs["is_grn"] == True

        conv_stf = False
        if g_sta[syn_tag].attrs["is_grn"]:
            conv_stf = True
            assert g_sta[syn_tag].attrs["origin_time"] == event_t0
            assert g_sta[dsyn_tag].attrs["origin_time"] == event_t0

        # time samples
        data_tb = obs.attrs["starttime"]
        data_fs = obs.attrs["sampling_rate"]
        data_dt = 1.0 / data_fs
        data_nt = obs.attrs["npts"]
        data_times = np.arange(data_nt, dtype=float) * data_dt

        for win in windows:
            # print(win["id"])
            # skip bad windows
            if not win["valid"]:
                msg = f"Window not measured for adj (WIN: {win}), SKIP"
                warnings.warn(msg)
                continue

            # window weight
            weight = win["weight"]

            # # skip window with zero weight
            # if np.isclose(weight, 0.0):
            #     msg = f"Window has near zero weight ({win['weight']}) (WIN: {win}), SKIP"
            #     warnings.warn(msg)
            #     continue

            weight_sum += weight

            # filter
            butter_N = win["butter_N"]
            butter_Wn = win["butter_Wn"]
            # fft
            npad = int(2 * data_fs / min(butter_Wn))
            nfft = scipy.fft.next_fast_len(data_nt + npad)
            freqs = np.fft.rfftfreq(nfft, d=data_dt)
            # window filter
            sos = scipy.signal.butter(
                butter_N, butter_Wn, "bandpass", fs=data_fs, output="sos"
            )
            _, filter_h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
            Fw = abs(filter_h)  # zero-phase filter response
            # time derivative in frequency domain
            Ft = 2j * np.pi * freqs
            # time window
            win_tb = UTCDateTime(win["starttime"])
            win_te = UTCDateTime(win["endtime"])
            win_taper = win["taper"]
            # window function
            b = win_tb - data_tb
            e = win_te - data_tb
            win_len = e - b
            width = win_len * min(win_taper, 0.4)
            win_c = [b, b + width, e - width, e]
            win_func = cosine_sac_taper(data_times, win_c)
            # polarity projection
            proj_matrix = win["proj_matrix"]
            #
            obs_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ, nfft), nfft)[:, :data_nt]
            obs_filt_win = np.dot(proj_matrix, obs_filt) * win_func
            norm_obs = np.sqrt(np.sum(obs_filt_win**2))

            f_syn = np.fft.rfft(syn_ENZ, nfft)
            if syn_to_vel:
                f_syn *= Ft
            f_dsyn_dm = {}
            for m, dsyn in dsyn_dm.items():
                f_dsyn = np.fft.rfft(dsyn, nfft)
                if syn_to_vel:
                    f_dsyn *= Ft
                f_dsyn_dm[m] = f_dsyn

            # search each model grid
            for istep in range(nstep):
                f_syn1 = np.zeros_like(f_syn)
                f_syn1 += f_syn
                for m, f_dsyn in f_dsyn_dm.items():
                    f_syn1 += dm[m][istep] * f_dsyn

                # perturbed source time function
                dt0 = 0.0
                if "dt0" in dm:
                    dt0 = dm["dt0"][istep]
                    f_tshift = np.exp(-2.0j * np.pi * freqs * dt0)
                    f_syn1 *= f_tshift

                dtau = 0.0
                if "dtau" in dm:
                    dtau = dm["dtau"][istep]
                if conv_stf:
                    f_src = stf_gauss_spectrum(freqs, event["tau"] + dtau)
                    f_syn1 *= f_src

                syn1_filt = np.fft.irfft(Fw * f_syn1, nfft)[:, :data_nt]
                syn1_filt_win = np.dot(proj_matrix, syn1_filt) * win_func
                norm_syn1 = np.sqrt(np.sum(syn1_filt_win**2))

                # normalized cc between obs and perturbed syn
                Nw = norm_obs * norm_syn1
                if Nw == 0:
                    msg = f"{net}.{sta}: obs data is zero for window({win})"
                    warnings.warn(msg)
                    continue
                cc = np.sum(obs_filt_win * syn1_filt_win) / Nw

                # weighted cc
                wcc_sum[istep] += weight * cc

    h5f.close()

    return wcc_sum, weight_sum


class Misfit(object):
    """Class managing all misfit windows

        Unit: kg,m,s
        Coordinate: ECEF Cartesian

    hdf5 data structure:
    /
    |- attrs: iteration_no., stage_type, inversion_type
    |         config: dictionary of parameters in misfit_par.yaml
    |         solver: dictionary of parameters in SEM/DATA/Par_file
    |
    |- source/  # point-source CMT representation,
    |   |- attrs: header, id,
    |   |         t0, tau, xs[3], mt[3,3]
    |   |         grad_xs[3], grad_mt[3,3]
    |   |         dCMT: {'dcmt1':{'dxs':dxs[3], 'dmt':dmt[3,3]}, 'dcmt2':{}}
    |   |         adj_strain[3,3,nt] # adjoint wavefield at source location
    |
    |- structure/
    |   |- attrs: perturbation_tags[:] (e.g. ['dvp', 'dvs'] )
    |
    |- solver/
    |   |- attrs: name, starttime, dt, nt
    |
    |- window (table)
    |   |- net, sta, name, time_win, filter_band, component, proj_matrix,
    |   |  cc0, AR0, cc_max_tshift, cc_max, ARmax,
    |   |  Asignal, Anoise, Asyn, SNR, weight,
    |
    |- waveforms/
    |   |- NET_STA/
    |   |   |- attrs: net,sta,loc,stla,stlo,stel,stdp,channels
    |   |   |         gcarc, az, baz, ttp
    |   |   |- DATA[ncha, nt]
    |   |   |   |- attrs: channels=[(code,az,dip), ...], filter
    |   |   |- SYN[ENU, nt]
    |   |   |- dSYN/
    |   |   |   |- dmodel[enu, nt]  # structure inversion
    |   |   |   |- dcmt1[enu,nt]   # source inversion
    |
    |- grid_search/
    |   |- source/
    |   |   |- dt0[nstep]
    |   |   |- dtau[nstep]
    |   |   |- dxs[nstep]
    |   |   |- dmt[nstep]
    |   |   |- wcc_sum[nstep]
    |   |   |   |- attrs: weight_sum
    |   |- structure/
    |   |   |- dvp[nstep] # step length relative to external model perturbation GLL files
    |   |   |- dvs[nstep]
    |   |   |- wcc_sum[nstep]
    |   |   |   |- attrs: weight_sum

    Methods:
      setup_event
      setup_station
      read_obs_syn
      read_perturbed_waveform
      delete_window
      add_window_body_wave
      add_window_surface_wave
      measure_adj
      window_quality_control (determine bad/OK, and window weight)
      output_adj
      make_cmt_dxs/dmt
      waveform_der_dxs/dmt
      waveform_der_dmodel
      cc_perturbed_seisomgram
      grid_cc

    NOTE:
      0. 1D Earth model: ak135

    """

    def __init__(self, h5_path):
        if os.path.isfile(h5_path):
            if pt.is_hdf5_file(h5_path):
                self.h5_path = h5_path
            else:
                raise FileExistsError(f"{h5_path} is not a valid HDF5 file!")
        else:
            with pt.open_file(h5_path, "w") as h5f:
                pass
        self.h5_path = h5_path

    def read_config_file(self, config_yaml):
        with open(config_yaml, "r") as file:
            self.config = yaml.safe_load(file)

        with pt.open_file(self.h5_path, "r+") as h5f:
            if "config" in h5f.root._v_attrs:
                msg = f"config field exists, overwrite!"
                warnings.warn(msg)
            h5f.root._v_attrs["config"] = self.config

    def read_solver_parfile(self, solver_parfile):
        """
        parfile: key = value
        """
        self.solver_param = pd.read_csv(
            solver_parfile,
            delimiter=r"\s*=\s*",
            header=None,
            comment="#",
            names=["key", "value"],
            dtype=dict(key=object, value=object),
            index_col=["key"],
        ).to_dict()["value"]

        with pt.open_file(self.h5_path, "r+") as h5f:
            if "solver" in h5f.root._v_attrs:
                msg = f"solver field exists, overwrite!"
                warnings.warn(msg)
            h5f.root._v_attrs["solver"] = self.solver_param

    def read_cmtsolution(self, cmt_file, ECEF=False):
        """cmt_file (str): CMTSOLUTION format file"""

        with open(cmt_file, "r") as f:
            lines = [x for x in f.readlines() if not (x.startswith("#"))]

        header = lines[0].split()
        year = header[1]
        month = header[2]
        day = header[3]
        hour = header[4]
        minute = header[5]
        second = header[6]

        lines = [x.split(":") for x in lines]
        event_id = lines[1][1].strip()
        time_shift = float(lines[2][1])

        # initialize pyproj objects
        # GPS_ELLPS = h5f.root._v_attrs["config"]["gps_ellps"]
        # ecef = pyproj.Proj(proj="geocent", ellps=GPS_ELLPS)
        # lla = pyproj.Proj(proj="latlong", ellps=GPS_ELLPS)
        gps2ecef = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")  # GPS to ECEF
        ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")  # ECEF to GPS

        if ECEF:
            tau = float(lines[3][1])
            x = float(lines[4][1])
            y = float(lines[5][1])
            z = float(lines[6][1])
            # convert from ECEF(meters) to lla
            lat, lon, alt = ecef2gps.transform(x, y, z)
            dep = -alt / 1000.0
        else:
            tau = float(lines[3][1]) / 1.628  # mimic triangle with gaussian
            lat = float(lines[4][1])
            lon = float(lines[5][1])
            dep = float(lines[6][1])
            # convert from lla to ECEF(meters)
            alt = -1000.0 * dep  # NOTE ignore local topography
            x, y, z = gps2ecef.transform(lat, lon, alt)
            # x, y, z = pyproj.transform(lla, ecef, lon, lat, alt)

        # centroid time: t0
        isotime = "{:s}-{:s}-{:s}T{:s}:{:s}:{:s}Z".format(
            year, month, day, hour, minute, second
        )
        t0 = UTCDateTime(isotime) + time_shift
        # modify origin time in header line to the centroid origin time
        header[1] = "{:04d}".format(t0.year)
        header[2] = "{:02d}".format(t0.month)
        header[3] = "{:02d}".format(t0.day)
        header[4] = "{:02d}".format(t0.hour)
        header[5] = "{:02d}".format(t0.minute)
        minisecond = round(t0.microsecond * 1e-3)
        header[6] = "{:02d}.{:03d}".format(t0.second, minisecond)
        header_line = " ".join(header)
        t0 = t0.replace(microsecond=minisecond * 1000)

        # moment tensor
        # ECEF=false: 1,2,3 -> r,theta,phi
        # ECEF=true:  1,2,3 -> x,y,z
        m11 = float(lines[7][1])
        m22 = float(lines[8][1])
        m33 = float(lines[9][1])
        m12 = float(lines[10][1])
        m13 = float(lines[11][1])
        m23 = float(lines[12][1])
        mt = np.array([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]])
        # transform from spherical to cartesian coordinate
        r = (x**2 + y**2 + z**2) ** 0.5
        theta = np.arccos(z / r)
        phi = np.arctan2(y, x)
        # rotation matrix
        sthe = np.sin(theta)
        cthe = np.cos(theta)
        sphi = np.sin(phi)
        cphi = np.cos(phi)
        # basis transform matrix: e_x,y,z = a * e_r,t,p
        mt_xyz = np.zeros((3, 3))
        mt_rtp = np.zeros((3, 3))
        a = np.array(
            [
                [sthe * cphi, cthe * cphi, -1.0 * sphi],
                [sthe * sphi, cthe * sphi, cphi],
                [cthe, -1.0 * sthe, 0.0],
            ]
        )
        if ECEF:
            mt_xyz = mt
            mt_rtp = np.dot(np.dot(np.transpose(a), mt), a)
        else:  # spherical coordinate
            # harvard cmt use dyn*cm, change to N*m
            mt_rtp = mt * 1.0e-7
            mt_xyz = np.dot(np.dot(a, mt_rtp), np.transpose(a))

        with pt.open_file(self.h5_path, "r+") as h5f:
            src_path = f"/source"
            if src_path in h5f:
                h5f.remove_node(src_path)
            tbl_src = h5f.create_table("/", "source", Source)

            src_row = tbl_src.row
            src_row["header"] = header_line
            src_row["id"] = event_id
            src_row["t0"] = t0.isoformat()
            src_row["tau"] = tau
            src_row["xs"] = [x, y, z]
            src_row["mt"] = mt_xyz
            src_row["mt_rtp"] = mt_rtp
            src_row["latitude"] = lat
            src_row["longitude"] = lon
            src_row["depth"] = dep  # km
            src_row.append()
            tbl_src.flush()

        # g_src._v_attrs["id"] = event_id
        # g_src._v_attrs["header"] = " ".join(header)
        # g_src._v_attrs["longitude"] = lon
        # g_src._v_attrs["latitude"] = lat
        # g_src._v_attrs["depth"] = dep
        # g_src._v_attrs["t0"] = t0
        # g_src._v_attrs["tau"] = tau
        # g_src._v_attrs["xs"] = np.array([x, y, z])
        # g_src._v_attrs["mt"] = mt_xyz
        # g_src._v_attrs["mt_rtp"] = mt_rtp

    def write_cmtsolution(self, out_file, ECEF=True):
        """
        Write out CMTSOLUTION file
        """
        with pt.open_file(self.h5_path, "r") as h5f:
            src_path = f"/source"
            if src_path not in h5f:
                msg = "no source information"
                raise ValueError(msg)
            else:
                tbl_src = h5f.get_node(src_path)
                if tbl_src.nrows == 0:
                    msg = "no source information"
                    raise ValueError(msg)

            event = tbl_src[0]
            event_id = event["id"].decode()
            header = event["header"].decode()
            tau = event["tau"]

        if ECEF:
            xs = event["xs"]
            mt = event["mt"]
            with open(out_file, "w") as fp:
                fp.write("%s\n" % (header))
                fp.write("%-18s %s\n" % ("event name:", event_id))
                fp.write("%-18s %+15.8E\n" % ("t0(s):", 0.0))
                fp.write("%-18s %+15.8E\n" % ("tau(s):", tau))
                fp.write("%-18s %+.3f\n" % ("x(m):", xs[0]))
                fp.write("%-18s %+.3f\n" % ("y(m):", xs[1]))
                fp.write("%-18s %+.3f\n" % ("z(m):", xs[2]))
                fp.write("%-18s %+15.8E\n" % ("Mxx(N*m):", mt[0, 0]))
                fp.write("%-18s %+15.8E\n" % ("Myy(N*m):", mt[1, 1]))
                fp.write("%-18s %+15.8E\n" % ("Mzz(N*m):", mt[2, 2]))
                fp.write("%-18s %+15.8E\n" % ("Mxy(N*m):", mt[0, 1]))
                fp.write("%-18s %+15.8E\n" % ("Mxz(N*m):", mt[0, 2]))
                fp.write("%-18s %+15.8E\n" % ("Myz(N*m):", mt[1, 2]))
        else:
            mt_rtp = event["mt_rtp"]
            with open(out_file, "w") as fp:
                fp.write("%s\n" % (header))
                fp.write("%-18s %s\n" % ("event name:", event_id))
                fp.write("%-18s %+15.8E\n" % ("time shift:", 0.0))
                fp.write("%-18s %+15.8E\n" % ("half duration:", tau * 1.628))
                fp.write("%-18s %+15.8E\n" % ("latitude:", event["latitude"]))
                fp.write("%-18s %+15.8E\n" % ("longitude:", event["longitude"]))
                fp.write("%-18s %+15.8E\n" % ("depth:", event["depth"]))
                fp.write("%-18s %+15.8E\n" % ("Mrr(dyn*cm):", mt_rtp[0, 0] * 1.0e7))
                fp.write("%-18s %+15.8E\n" % ("Mtt(dyn*cm):", mt_rtp[1, 1] * 1.0e7))
                fp.write("%-18s %+15.8E\n" % ("Mpp(dyn*cm):", mt_rtp[2, 2] * 1.0e7))
                fp.write("%-18s %+15.8E\n" % ("Mrt(dyn*cm):", mt_rtp[0, 1] * 1.0e7))
                fp.write("%-18s %+15.8E\n" % ("Mrp(dyn*cm):", mt_rtp[0, 2] * 1.0e7))
                fp.write("%-18s %+15.8E\n" % ("Mtp(dyn*cm):", mt_rtp[1, 2] * 1.0e7))

    def read_channel_file(self, channel_file, station_file=None):
        """
        channel_file : str
            FDSN-station text format file at channel level

        station (table)
            net,sta,loc,bi,components,stla,stlo,stel,stdp
        """
        # config = h5f.root._v_attrs["config"]

        # if "/source" not in h5f:
        #     raise KeyError("/source does not exist, run read_cmtsolutin first")
        # tbl_src = h5f.get_node("/source")
        # if tbl_src.nrows == 0:
        #     msg = "no source information"
        #     raise Exception(msg)
        # event = tbl_src[-1]
        # origin_time = UTCDateTime(event["t0"])

        h5f = pt.open_file(self.h5_path, "r+")
        if "/channel" in h5f:
            h5f.remove_node("/", "channel")
        tbl_chan = h5f.create_table("/", "channel", Channel)
        chan_row = tbl_chan.row

        select_station = False
        if station_file:
            try:
                stations_df = pd.read_csv(
                    station_file, sep=r"\s+", header=None, comment="#"
                )
                select_station = True
            except:
                msg = f"failed to read station file, ignored ({station_file})"
                warnings.warn(msg)
                select_station = False

        with open(channel_file, "r") as f:
            lines = [
                x.replace("\n", "").split("|")
                for x in f.readlines()
                if not (x.startswith("#"))
            ]

        # sort channels based on station_id (net,sta,loc,bi) bi=band,instrument
        for l in lines:
            l = [x.strip() for x in l]

            net, sta, loc, cha = l[0], l[1], l[2], l[3]
            stla, stlo, stel, stdp, az, dip = (float(x) for x in l[4:10])

            if select_station:
                match = (stations_df.iloc[:, 0] == net) & (
                    stations_df.iloc[:, 1] == f"{sta}.{loc}"
                )
                if not any(match):
                    # msg = f"no match in station file ({station_file}), skip ({l})"
                    continue

            if len(cha) != 3:
                msg = f"wrong length of channel name (!=3) in channel info, skip ({l})"
                warnings.warn(msg)
                continue

            if l[15]:
                tb = UTCDateTime(l[15])
            else:
                msg = f"no starttime in channel info, skip ({l})"
                warnings.warn(msg)
                continue
            if l[16]:
                te = UTCDateTime(l[16])
            else:
                # msg = f"empty endtime in ({l})"
                # warnings.warn(msg)
                te = None  # UTCDateTime(empty_endtime)
            # date1 = [int(a) for a in re.sub("\D", " ", x[15]).split()]
            # date2 = [int(a) for a in re.sub("\D", " ", x[16]).split()]
            # t1 = (
            #     UTCDateTime(date1[0], date1[1], date1[2])
            #     + 60.0 * (60.0 * date1[3] + date1[4])
            #     + date1[5]
            # )
            # t2 = (
            #     UTCDateTime(date2[0], date2[1], date2[2])
            #     + 60.0 * (60.0 * date2[3] + date2[4])
            #     + date2[5]
            # )
            chan_row["network"] = net
            chan_row["station"] = sta
            chan_row["location"] = loc
            chan_row["channel"] = cha
            chan_row["latitude"] = stla
            chan_row["longitude"] = stlo
            chan_row["elevation"] = stel
            chan_row["depth"] = stdp
            chan_row["azimuth"] = az
            chan_row["dip"] = dip
            chan_row["starttime"] = tb.isoformat()
            chan_row["endtime"] = te.isoformat() if te else ""

            chan_row.append()

        tbl_chan.flush()

        h5f.close()

    def read_data_h5(self, h5_file, station_file=None):
        """
        Read in observed seismograms from hdf5 file.
        """
        h5f = pt.open_file(self.h5_path, "r+")

        config = h5f.root._v_attrs["config"]
        time_before_first_arrival = config["data"]["time_before_first_arrival"]
        if not isinstance(time_before_first_arrival, list):
            time_before_first_arrival = [
                time_before_first_arrival,
            ]
        # taper_width = config['data']['taper_width']
        time_after_origin = config["data"]["time_after_origin"]
        if not isinstance(time_after_origin, list):
            time_after_origin = [
                time_after_origin,
            ]
        obs_tag = config["data"]["tag"]

        if "/source" not in h5f:
            raise KeyError("/source does not exist, run read_cmtsolutin first")
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])

        geod = pyproj.Geod(ellps=config["gps_ellps"])
        R_earth = (geod.a + geod.b) / 2
        taup_model = TauPyModel(model=config["taup_model"])

        if "/waveform" not in h5f:
            g_wav = h5f.create_group("/", "waveform")
        else:
            g_wav = h5f.get_node("/waveform")

        if "/channel" not in h5f:
            raise KeyError("/channel does not exist, run read_channel_file first")
        tbl_chan = h5f.get_node("/channel")

        select_station = False
        if station_file:
            try:
                stations_df = pd.read_csv(
                    station_file, sep=r"\s+", header=None, comment="#"
                )
                select_station = True
            except:
                msg = f"failed to read station file ({station_file})"
                warnings.warn(msg)

        # for storing data, /waveforms/NET_STA/DATA_DISP[nchan,nt]
        h5_atom = pt.Atom.from_dtype(np.dtype(np.float32))
        h5_filters = pt.Filters(complevel=3, complib="zlib")

        with pt.open_file(h5_file, "r") as obs_h5f:
            for g_sta_obs in obs_h5f.root:
                # print(g_sta_obs._v_name)

                sta_attrs = g_sta_obs._v_attrs
                net = sta_attrs["network"]
                sta = sta_attrs["station"]
                loc = sta_attrs["location"]

                if obs_tag not in g_sta_obs:
                    msg = f"{g_sta_obs._v_name}: {obs_tag} does not exist, skip"
                    warnings.warn(msg)
                    continue
                data = g_sta_obs[obs_tag]
                data_fs = data.attrs["sampling_rate"]
                data_dt = 1.0 / data_fs
                data_starttime = data.attrs["starttime"]
                data_npts = data.attrs["npts"]
                data_endtime = data_starttime + (data_npts - 1) * data_dt
                data_channels = data.attrs["channels"]

                if select_station:
                    match = (stations_df.iloc[:, 0] == net) & (
                        stations_df.iloc[:, 1] == f"{sta}.{loc}"
                    )
                    if not any(match):
                        # msg = f"not in station select list, skip ({net}.{sta})"
                        # warnings.warn(msg)
                        continue

                channels = tbl_chan.read_where(
                    f'(network == b"{net}") & (station == b"{sta}") & (location == b"{loc}")'
                )
                channels = [
                    cha
                    for cha in channels
                    if (UTCDateTime(cha["starttime"]) < data_starttime)
                    and (
                        (not cha["endtime"])
                        or (UTCDateTime(cha["endtime"]) > data_endtime)
                    )
                ]

                if not channels:
                    msg = f"no channel info for {net}.{sta}.{loc}, ignored"
                    warnings.warn(msg)
                    continue

                # channel_names = [cha["channel"] for cha in channels]
                # match = [cha["name"] in channel_names for cha in data_channels]
                # if not all(match):
                #     msg = f"missing data channel in {channels}, ignored"
                #     warnings.warn(msg)
                #     continue

                # update az/dip from channel_file
                for cha in data_channels:
                    cha1 = [c for c in channels if c["channel"] == cha["name"]]
                    if len(cha1) == 0:
                        msg = f"no channel info found for {cha["name"]} ({g_sta_obs._v_name}),  use az/dip given in the data"
                        warnings.warn(msg)
                        continue
                    if len(cha1) > 1:
                        msg = f"multiple channel info found ({cha1}),  use the first returned values"
                        warnings.warn(msg)
                    cha["azimuth"] = cha1[0]["azimuth"]
                    cha["dip"] = cha1[0]["dip"]

                coords = [
                    (cha["latitude"], cha["longitude"], cha["elevation"], cha["depth"])
                    for cha in channels
                ]
                if len(set(coords)) != 1:
                    msg = f"inconsistent coordinates in ({channels}), ignored"
                    warnings.warn(msg)
                    continue
                stla, stlo, stel, stdp = coords[0]

                # check vertical component
                Z_comp = [
                    (cha["name"], cha["azimuth"], cha["dip"])
                    for cha in data_channels
                    if cha["name"].decode()[-1] == "Z"
                ]
                if len(Z_comp) != 1 or abs(Z_comp[0][2]) != 90.0:
                    msg = f"problematic vertical channel, {Z_comp}"
                    warnings.warn(msg)
                    continue

                # check horizontal components
                H_comp = [
                    (cha["name"], cha["azimuth"], cha["dip"])
                    for cha in data_channels
                    if cha["name"].decode()[-1] != "Z"
                ]
                if len(H_comp) == 0:
                    msg = f"no horizontal channels found for {g_sta_obs._v_name}"
                    warnings.warn(msg)
                elif (
                    len(H_comp) != 2
                    or abs(H_comp[0][2]) != 0.0
                    or abs(H_comp[1][2]) != 0.0
                    or abs(np.cos(np.deg2rad(H_comp[0][1] - H_comp[1][1]))) > 0.1
                ):
                    msg = f"problematic horizontal channels for {g_sta_obs._v_name}: {H_comp}, ignored"
                    warnings.warn(msg)
                    H_comp = []

                # calculate first arrival time
                az, baz, dist = geod.inv(
                    event["longitude"], event["latitude"], stlo, stla
                )
                # az, baz in range [0, 360)
                az = az % 360
                baz = baz % 360
                dist_degree = np.rad2deg(dist / R_earth)
                evdp_km = event["depth"]
                if evdp_km < 0.0:
                    evdp_km = 0.0
                ttp = taup_model.get_travel_times(
                    source_depth_in_km=evdp_km,
                    distance_in_degree=dist_degree,
                    phase_list=["ttp"],
                )
                first_arrtime = event_t0 + min([arr.time for arr in ttp])

                # check if data time range includes the specified minimum time window
                t0 = first_arrtime - min(time_before_first_arrival)
                t1 = event_t0 + min(time_after_origin)
                if data_starttime > t0 or data_endtime < t1:
                    msg = f"{g_sta_obs._v_name}: timespan [{data_starttime}, {data_endtime}] does not cover required [{t0}, {t1}], skip"
                    warnings.warn(msg)
                    continue

                # create station group, e.g. /waveforms/NET_STA
                station_name = f"{net}_{sta}"
                if station_name in g_wav:
                    warnings.warn(f"{station_name} alreay exists, overwrite.")
                    g_sta = h5f.get_node(g_wav, station_name)
                else:
                    g_sta = h5f.create_group(g_wav, station_name)
                for attr_name in g_sta_obs._v_attrs._f_list("user"):
                    g_sta._v_attrs[attr_name] = g_sta_obs._v_attrs[attr_name]
                g_sta._v_attrs["first_arrtime"] = first_arrtime
                g_sta._v_attrs["azimuth"] = az
                g_sta._v_attrs["back_azimuth"] = baz
                g_sta._v_attrs["dist_degree"] = dist_degree

                # link waveform data, e.g. /waveforms/NET_STA/DATA_DISP --> data.h5:/NET_STA/DATA_DISP
                if obs_tag in g_sta:
                    msg = f"{obs_tag} alreay exists, overwrite."
                    warnings.warn(msg)
                    h5f.remove_node(g_sta, obs_tag)
                # h5f.create_external_link(g_sta, obs_tag, data)

                # cut data /waveforms/NET_STA/DATA_DISP
                t0 = first_arrtime - max(time_before_first_arrival)
                t1 = event_t0 + max(time_after_origin)
                t0 = max(t0, data_starttime)  # cut window begins no earlier than t0
                t1 = min(t1, data_endtime)  # cut window ends no later than t1
                idx0 = int((t0 - data_starttime) * data_fs)
                idx0 = max(idx0, 0)
                idx1 = int((t1 - data_starttime) * data_fs)
                idx1 = min(idx1, data_npts)
                npts = idx1 - idx0
                starttime = data_starttime + idx0 * data_dt
                nchan = data.shape[0]
                shape = (nchan, npts)
                ca = h5f.create_carray(
                    g_sta, obs_tag, h5_atom, shape, filters=h5_filters
                )
                # times = np.arange(npts) * data_dt
                # e = (npts - 1) * data_dt
                # win_c = [0, taper_width, e -taper_width, e]
                # win_func = cosine_sac_taper(times, win_c)
                # data_win = data[:, idx0:idx1] # * win_func
                # ca[:] = (data_win - np.mean(data_win, axis=-1, keepdims=True)) * win_func
                ca[:] = data[:, idx0:idx1]
                ca.attrs["starttime"] = starttime
                ca.attrs["sampling_rate"] = data_fs
                ca.attrs["npts"] = npts
                ca.attrs["channels"] = data_channels
                ca.attrs["filter"] = data.attrs["filter"]
                # ca.attrs['taper'] = win_c
                ca.attrs["type"] = data.attrs["type"]

                # if "md5" not in g_sta._v_attrs:
                #     g_sta._v_attrs["md5"] = {}
                # g_sta._v_attrs["md5"][obs_tag] = hashlib.md5(data).hexdigest()

                h5f.flush()

        h5f.close()

    def read_syn_sac(
        self,
        sac_dir,
        is_grn=False,
        is_diff=False,
        tag_diff="diff",
        # diff_src="/source/dxs_dmt",
    ):
        """
        sac file naming rule: NET.STA.LOC.CHA , where CHA consists of band_code(e.g. BH or MX) + orientation[E|N|Z]
        syn should be sac files and have sac header b and o set correctly.
        /waveforms/NET_STA/[SYN|GRN]_DISP[3, nt]
        """
        syn_components = ["E", "N", "Z"]

        h5f = pt.open_file(self.h5_path, "r+")

        config = h5f.root._v_attrs["config"]

        # time_before_first_arrival_as_noise = config["data"][
        #     "time_before_first_arrival_as_noise"
        # ]
        # time_before_first_arrival = config["data"]["time_before_first_arrival"]

        data_tag = config["data"]["tag"]
        if is_diff:
            # syn_tag = config["syn"]["tag_diff"]
            dsyn_grp = config["syn"]["dsyn_grp"]
            dsyn_tag = tag_diff
            syn_tag = config["syn"]["tag"]
        else:
            syn_tag = config["syn"]["tag"]
        syn_type = config["syn"]["type"]
        # left_pad = config["syn"]["left_pad"]
        # right_pad = config["syn"]["right_pad"]
        syn_band_code = config["syn"]["band_code"]
        syn_sac_suffix = config["syn"]["sac_suffix"]
        # syn_is_grn = config["syn"]["is_grn"]

        # if left_pad < 0:
        #     warnings.warn("[arg] left_pad < 0, set to 0.0")
        #     left_pad = 0
        # if right_pad < 0:
        #     warnings.warn("[arg] right_pad < 0, set to 0.0")
        #     right_pad = 0
        # if time_before_first_arrival_as_noise < 0:
        #     time_before_first_arrival_as_noise = 50.0
        #     warnings.warn("[arg] time_before_first_arrival_as_noise < 0, set to 50")

        if "/source" not in h5f:
            raise KeyError("/source does not exist, run read_cmtsolutin first")
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])

        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        for g_sta in g_wav:
            net = g_sta._v_attrs["network"]
            sta = g_sta._v_attrs["station"]
            loc = g_sta._v_attrs["location"]
            # stla = g_sta.attrs["latitude"]
            # stlo = g_sta.attrs["longitude"]
            # stel = g_sta.attrs["elevation"]
            # stdp = g_sta.attrs["depth"]
            # first_arrtime = g_sta.attrs["first_arrtime"]

            if data_tag in g_sta:
                obs_data = g_sta[data_tag]
            else:
                msg = f"{data_tag} does not exist under {g_sta._v_file.filename}:{g_sta._v_pathname}, skip"
                warnings.warn(msg)
                continue

            if is_diff:
                if syn_tag in g_sta:
                    syn = g_sta[syn_tag]
                    if is_grn:
                        if not syn.attrs["is_grn"]:
                            msg = f"{syn._v_pathname}.attrs[is_grn]==False but is_grn is set to True!"
                            raise ValueError(msg)
                        assert syn.attrs["origin_time"] == event_t0
                else:
                    msg = f"{syn_tag} does not exist under {g_sta._v_pathname}, skip"
                    warnings.warn(msg)
                    continue

            # try:
            #     obs_data = obs(mode="r")
            # except Exception as e:
            #     msg = f"failed to dereference link {obs.target} (ERROR: {e})"
            #     warnings.warn(msg)
            #     continue
            obs_fs = obs_data.attrs["sampling_rate"]
            obs_dt = 1.0 / obs_fs
            obs_tb = obs_data.attrs["starttime"]
            obs_nt = obs_data.attrs["npts"]
            obs_te = obs_tb + (obs_nt - 1) * obs_dt
            obs_filters = obs_data.attrs["filter"]

            # butter_N = obs_filter["N"]
            # butter_Wn = obs_filter["Wn"]
            # obs_taper = obs_data.attrs['taper']
            # times = np.arange(obs_nt) * obs_dt
            # e = (obs_nt - 1) * obs_dt
            # win_func = cosine_sac_taper(times, obs_taper)

            # read in synthetic seismograms

            syn_st = Stream()
            syn_sac_files = [
                os.path.join(
                    sac_dir, f"{net}.{sta}.{loc}.{syn_band_code}{comp}{syn_sac_suffix}"
                )
                for comp in syn_components
            ]
            read_sac_ok = True
            for sac_file in syn_sac_files:
                if not os.path.isfile(sac_file):
                    msg = f"{sac_file} does not exist"
                    warnings.warn(msg)
                    read_sac_ok = False
                    break
                try:
                    st1 = read(sac_file)
                except Exception as e:
                    msg = f"Error reading {sac_file} (Details: {e})"
                    warnings.warn(msg)
                    read_sac_ok = False
                    break
                if len(st1) != 1:
                    msg = f"more than one trace in {sac_file} ({st1})"
                    warnings.warn(msg)
                    read_sac_ok = False
                    break
                syn_st += st1
            if not read_sac_ok:
                continue

            if len(syn_st) != 3:
                msg = f"Not exactly 3 components found ! ({syn_st}), skip"
                warnings.warn(msg)
                continue
            if not is_equal(
                [(tr.stats.starttime, tr.stats.delta, tr.stats.npts) for tr in syn_st]
            ):
                msg = f"{syn_sac_files} do not have equal time samples."
                raise Exception(msg)

            tr = syn_st[0]
            solver_dt = tr.stats.delta
            solver_fs = 1.0 / solver_dt
            solver_nt = tr.stats.npts
            solver_tb = tr.stats.starttime
            solver_te = solver_tb + (solver_nt - 1) * solver_dt

            assert solver_fs > obs_fs

            # commented out to enforce consistency between sac o time and event t0 even for Green's function
            # if is_grn:  # for green function, use event['t0'] for the o time
            #    solver_tb = event_t0 - (tr.stats.sac["o"] - tr.stats.sac["b"])

            # check if the origin time in sac files is consistent with event['t0']
            syn_origin_time = solver_tb + (tr.stats.sac["o"] - tr.stats.sac["b"])
            if abs(syn_origin_time - event_t0) > 1.0e-6:
                err = f"{g_sta._v_name}: Inconsistent origin time between sac headers ({syn_origin_time}) and event['t0'] ({event_t0})"
                raise Exception(err)

            if is_diff:
                if (
                    solver_tb != syn.attrs["solver_tb"]
                    or solver_nt != syn.attrs["solver_nt"]
                    or solver_dt != syn.attrs["solver_dt"]
                ):
                    msg = f"perturbed synthetics ({syn_sac_files}) do not have equal time samples as the initial synthetics."
                    raise Exception(msg)

            # use_tb = max(first_arrtime - left_pad, solver_tb - left_pad, obs_tb)
            # use_te = min(solver_te, obs_te)

            # synthetic seismograms will be interpolated onto the same time samples as observed data
            # obs_idx0 = max(0, int((solver_tb - obs_tb) * obs_fs + 1))
            # obs_idx1 = min(obs_nt, int((solver_te - obs_tb) * obs_fs))

            nt_lpad = nt_rpad = 0
            if solver_tb > obs_tb:
                nt_lpad = int((solver_tb - obs_tb) * solver_fs + 5)
            if solver_te < obs_te:
                nt_rpad = int((obs_te - solver_te) * solver_fs + 5)

            # syn_fs = obs_fs
            # syn_dt = obs_dt
            # syn_tb = obs_tb + obs_idx0 * obs_dt
            # syn_nt = obs_idx1 - obs_idx0
            # for interpolation
            # ii = 0
            # if solver_tb > syn_tb:
            #     ii = int((solver_tb - syn_tb) * obs_fs + 1)
            # t0 = syn_tb + ii * obs_dt
            # nt = syn_nt - ii
            # t0 = obs_tb + obs_idx0 * obs_dt
            # nt = obs_idx1 - obs_idx0

            # apply the same filter used on the observed seismograms
            t0 = solver_tb - nt_lpad * solver_dt
            nt = solver_nt + nt_lpad + nt_rpad

            min_cutoff_freq = min([filter["Wn"] for filter in obs_filters])
            pad_length = max(20.0, 2.0 / min_cutoff_freq)
            npad = int(pad_length * solver_fs)
            nfft = scipy.fft.next_fast_len(nt + 2 * npad)  # pad both ends of the signal
            freqs = np.fft.rfftfreq(nfft, d=solver_dt)
            filter_h = np.ones_like(freqs)
            for filter in obs_filters:
                sos = scipy.signal.butter(
                    filter["N"],
                    filter["Wn"],
                    btype=filter["type"],
                    fs=solver_fs,
                    output="sos",
                )
                _, h = scipy.signal.freqz_sos(sos, worN=freqs, fs=solver_fs)
                filter_h *= abs(h)
            # avoid filter response from edge
            valid_te = solver_te - 1.0 / min_cutoff_freq

            #
            taper = cosine_taper(
                np.arange(nfft),
                [0, npad + nt_lpad, npad + nt_lpad + solver_nt, nfft - 1],
            )

            if _DEBUG:
                syn_st1 = syn_st.copy()

            for i in range(3):
                tr = syn_st[i]
                # st_tmp = Stream()
                # st_tmp += tr.copy()
                assert tr.id[-1] == syn_components[i]
                # x = np.zeros(nt)
                # x[nt_lpad : (solver_nt + nt_lpad)] = tr.data
                # extend synthetic data by edge values from both ends
                x = (
                    np.pad(
                        tr.data,
                        (npad + nt_lpad, nfft - npad - nt_lpad - solver_nt),
                        "edge",
                    )
                    * taper
                )
                x = np.fft.irfft(np.fft.rfft(x, nfft) * abs(filter_h))[
                    npad : (npad + nt)
                ]

                # # NOTE: different precision between UTCDateTime.timestamp and (t1 -t2)
                # #       UTCDateTime internally stores time in nanoseconds as an integer
                # #       t.timestamp = float64(t.ns * 1.e-9)
                # #       t1 - t2 = round((t1.ns - t2.ns)*1.e-9, UTCDateTime.precision)
                # #       to make things consistent, set UTCDateTime.DEFAULT_PRECISION = 9 after import
                # #       (t1.timestamp - t2.timestamp) used in tr.interpolate has less precision than (t1-t2)
                # tr.stats.starttime = t0
                # tr.data = x
                # # st_tmp += tr.copy()
                # tr.interpolate(
                #     obs_fs, starttime=obs_tb, npts=obs_nt, method="lanczos", a=20
                # )

                # call lanczos_interpolation directly
                tr.data = lanczos_interpolation(
                    x, 0, solver_dt, obs_tb - t0, obs_dt, obs_nt, a=20
                )
                tr.stats.starttime = obs_tb
                tr.stats.sampling_rate = obs_fs

                # st_tmp += tr.copy()
                # st_tmp[0].stats.location = '00'
                # st_tmp[1].stats.location = '01'
                # st_tmp[2].stats.location = '02'
                # st_tmp.plot()

            if _DEBUG:
                for tr in syn_st:
                    tr.stats.location = "SS"
                syn_st1.extend(syn_st)
                print("DEBUG:")
                print(syn_st1)
                syn_st1.plot()

            # store synthetic waveforms, e.g. /NET_STA/SYN_DISP[0:nchan, 0:npts]
            atom = pt.Atom.from_dtype(np.dtype(np.float32))
            filters = pt.Filters(complevel=3, complib="zlib")
            shape = (3, obs_nt)

            if is_diff:
                if dsyn_grp not in g_sta:
                    g_dsyn = h5f.create_group(g_sta, dsyn_grp)
                else:
                    g_dsyn = h5f.get_node(g_sta, dsyn_grp)
                if dsyn_tag in g_dsyn:
                    msg = f"{dsyn_tag} exists in {g_dsyn._v_pathname}, overwrite!"
                    warnings.warn(msg)
                    h5f.remove_node(g_dsyn, dsyn_tag)
                ca = h5f.create_carray(g_dsyn, dsyn_tag, atom, shape, filters=filters)
            else:
                if syn_tag in g_sta:
                    msg = f"{syn_tag} exists in {g_sta._v_pathname}, overwrite!"
                    warnings.warn(msg)
                    h5f.remove_node(g_sta, syn_tag)
                ca = h5f.create_carray(g_sta, syn_tag, atom, shape, filters=filters)

            ca[:] = 0
            for i in range(3):
                ca[i, :] = np.array(syn_st[i].data, dtype=np.float32)  # * win_func
            if is_diff:
                ca[:] -= syn[:]
            dtype = [("name", "S3"), ("azimuth", float), ("dip", float)]
            channels = [
                (syn_st[0].stats.channel, 90, 0),
                (syn_st[1].stats.channel, 0, 0),
                (syn_st[2].stats.channel, 0, -90),
            ]
            ca.attrs["channels"] = np.array(channels, dtype=dtype)
            ca.attrs["starttime"] = obs_tb
            ca.attrs["sampling_rate"] = obs_fs
            ca.attrs["npts"] = obs_nt
            ca.attrs["solver_tb"] = solver_tb
            ca.attrs["solver_dt"] = solver_dt
            ca.attrs["solver_nt"] = solver_nt
            ca.attrs["valid_te"] = valid_te
            ca.attrs["type"] = syn_type
            ca.attrs["is_grn"] = is_grn
            if is_grn:
                ca.attrs["origin_time"] = event_t0
        h5f.close()

    def setup_windows(self, window_yaml=None):
        """ """
        h5f = pt.open_file(self.h5_path, "r+")

        config = h5f.root._v_attrs["config"]

        window_cfg = None
        if window_yaml is not None:
            try:
                with open(window_yaml, "r") as file:
                    data = yaml.safe_load(file)
                window_cfg = data["window"]
            except Exception as e:
                msg = f"failed to read in {window_yaml} (ERROR: {e})"
                warnings.warn(msg)

        if window_cfg is not None:
            if "window" in config:
                msg = f"root._v_attrs[window] exists in {h5f.filename}, overwrite!"
                warnings.warn(msg)    
            config["window"] = window_cfg
        elif "window" in config:
            window_cfg = config["window"]
        else:
            msg = "failed to get window configuration"
            raise AssertionError(msg)

        win_taper = config["window_taper"]
        win_ids = set([w["id"] for w in window_cfg])
        if len(win_ids) != len(window_cfg):
            msg = "window's id is not unique in {window_cfg}"
            raise AssertionError(msg)

        cache = OrderedDict()
        taup_model = TauPyModel(model=config["taup_model"], cache=cache)

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])
        evdp_km = event["depth"]
        if evdp_km < 0.0:
            evdp_km = 0.0

        obs_tag = config["data"]["tag"]
        syn_tag = config["syn"]["tag"]

        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        if "/window" in h5f:
            tbl_win = h5f.remove_node("/", "window")
        tbl_win = h5f.create_table("/", "window", Window)
        win_row = tbl_win.row

        for g_sta in g_wav:
            # print(f"setup windows for {g_sta._v_name}")
            meta = g_sta._v_attrs
            net = meta["network"]
            sta = meta["station"]
            baz = meta["back_azimuth"]
            gcarc = meta["dist_degree"]

            if obs_tag not in g_sta:
                msg = f"{obs_tag} does not exist under {g_sta._v_file.filename}:{g_sta._v_pathname}, skip"
                warnings.warn(msg)
                continue
            if syn_tag not in g_sta:
                msg = f"{syn_tag} does not exist under {g_sta._v_file.filename}:{g_sta._v_pathname}, skip"
                warnings.warn(msg)
                continue

            obs = g_sta[obs_tag]
            syn = g_sta[syn_tag]
            solver_tb = syn.attrs["solver_tb"]
            solver_nt = syn.attrs["solver_nt"]
            solver_dt = syn.attrs["solver_dt"]
            solver_te = solver_tb + solver_nt * solver_dt
            solver_valid_te = syn.attrs["valid_te"]
            obs_tb = obs.attrs["starttime"]
            obs_nt = obs.attrs["npts"]
            obs_fs = obs.attrs["sampling_rate"]
            obs_te = obs_tb + obs_nt / obs_fs

            # valid_tb = max(solver_tb, obs_tb)
            valid_tb = obs_tb  # since synthetic is zero at any time before origin time
            # valid_te = min(solver_te, obs_te)
            valid_te = min(solver_valid_te, obs_te)

            phase_list = [
                phase.strip()
                for win in window_cfg
                if win["type"] == "body"
                for phase in win["phase"].split(",")
            ]
            phase_list = list(set(phase_list))
            all_arrivals = taup_model.get_travel_times(
                source_depth_in_km=evdp_km,
                distance_in_degree=gcarc,
                phase_list=phase_list,
            )
            p_arrivals = taup_model.get_travel_times(
                source_depth_in_km=evdp_km,
                distance_in_degree=gcarc,
                phase_list=["p", "P"],
            )
            min_ttp = min([arr.time for arr in p_arrivals])
            s_arrivals = taup_model.get_travel_times(
                source_depth_in_km=evdp_km,
                distance_in_degree=gcarc,
                phase_list=["s", "S"],
            )
            min_tts = min([arr.time for arr in s_arrivals])

            obs_channels = obs.attrs["channels"]
            obs_Zchan = [cha for cha in obs_channels if cha["name"].decode()[-1] == "Z"]
            obs_Hchan = [cha for cha in obs_channels if cha["name"].decode()[-1] != "Z"]
            if len(obs_Zchan) != 1 or obs_Zchan[0]["dip"] != -90:
                msg = (
                    f"{g_sta._v_name} has problematic Z channel info: {obs_Zchan}, skip"
                )
                raise AssertionError(msg)
            no_Hchan = False
            if len(obs_Hchan) == 0:
                msg = f"{g_sta._v_name} has no horizontal channels"
                warnings.warn(msg)
                no_Hchan = True
            else:
                assert len(obs_Hchan) == 2

            for win in window_cfg:
                win_id = win["id"]
                win_type = win["type"]
                win_phase = win["phase"]
                cmpnm = win["cmp"]
                win_bp = win["bp"] # filter parameter defined as [order, freq_min, freq_max]
                bp_long_period = 1.0 / min(win_bp[1:])
                evdp_limit = win["evdp"]
                gcarc_limit = win["gcarc"]
                win_weight = win["weight"]

                if no_Hchan and cmpnm != "Z":
                    msg = (
                        f"{g_sta._v_name}: no horizontal channels, skip window ({win})"
                    )
                    warnings.warn(msg)
                    continue

                if evdp_km < min(evdp_limit) or evdp_km > max(evdp_limit):
                    continue
                if gcarc < min(gcarc_limit) or gcarc > max(gcarc_limit):
                    continue

                if win_type == "body":
                    twin = win["twin"]
                    phase_list = [phase.strip() for phase in win_phase.split(",")]
                    arrivals = [arr for arr in all_arrivals if arr.name in phase_list]
                    # arrivals = taup_model.get_travel_times(
                    #     source_depth_in_km=evdp_km,
                    #     distance_in_degree=gcarc,
                    #     phase_list=phase_list,
                    # )
                    if not arrivals:
                        msg = f"{win_phase} not found for (evdp={evdp_km},gcarc={gcarc}), skip"
                        warnings.warn(msg)
                        continue
                    ttime = np.array([arr.time for arr in arrivals])
                    tb = min(ttime) + min(twin)
                    te = max(ttime) + max(twin)
                    # minlen = max(twin) - min(twin)
                elif win_type == "Pwin":
                    twin, dts = win["twin"], win["dts"]
                    tb = min_ttp + min(twin)
                    te = max(min_ttp + max(twin), min_tts + dts)
                elif win_type == "Swin":
                    twin, min_slowness = win["twin"], win["minslw"]
                    tb = min_tts + min(twin)
                    te = max(min_tts + max(twin), gcarc * min_slowness)
                elif win_type == "surf":
                    swin = win["swin"]
                    smin, smax, minlen = swin
                    tb = gcarc * smin
                    te = gcarc * smax
                    if (te - tb) < minlen:
                        tpad = 0.5 * (minlen - (te - tb))
                        tb -= tpad
                        te += tpad
                else:
                    msg = f"unknown window type: {win_type}, skip"
                    warnings.warn(msg)
                    continue

                # win_tb = max(valid_tb + bp_long_period, event_t0 + tb)
                # win_te = min(
                #     valid_te - bp_long_period, event_t0 + te
                # )  # +/-bp_long_period: to avoid filter edge effect in the later bandpass process
                win_tb = max(valid_tb, event_t0 + tb)
                win_te = min(valid_te, event_t0 + te)  

                # print(solver_te, obs_te, solver_valid_te, win_te)

                required_winlen = tb - te
                if (win_te - win_tb) < 0.5 * required_winlen:
                    msg = f"window ({win})'s valid length ({win_tb}-{win_te} is less than half of the required length ({required_winlen})"
                    warnings.warn(msg)
                    continue

                # window component
                if cmpnm == "Z":  # vertcal component
                    cmpaz = 0.0
                    cmpdip = -90.0
                elif cmpnm == "R":  # radial component
                    cmpaz = (baz + 180.0) % 360.0
                    cmpdip = 0.0
                elif cmpnm == "T":  # tangential component (TRZ: right-hand convention)
                    cmpaz = (baz - 90.0) % 360.0
                    cmpdip = 0.0
                elif cmpnm == "H":  # horizontal particle motion
                    cmpaz = float("nan")
                    cmpdip = 0.0
                elif cmpnm == "F":  # 3-d particle motion
                    cmpaz = float("nan")
                    cmpdip = float("nan")
                else:
                    warnings.warn("[WARN] %s: unrecognized component, SKIP." % (cmpnm))
                    continue

                #
                win_row["id"] = win_id
                win_row["network"] = net
                win_row["station"] = sta
                win_row["type"] = win_type
                win_row["phase"] = win_phase
                win_row["cmpnm"] = cmpnm
                win_row["cmpaz"] = cmpaz
                win_row["cmpdip"] = cmpdip
                # win_row["twin"] = twin
                win_row["butter_N"] = win_bp[0]
                win_row["butter_Wn"] = win_bp[1:]
                win_row["starttime"] = win_tb.isoformat()
                win_row["endtime"] = win_te.isoformat()
                win_row["taper"] = win_taper
                win_row["weight"] = win_weight
                win_row["valid"] = True

                win_row.append()

        tbl_win.flush()

        h5f.close()

    def __obsolete_get_obs_ENZ(self, g_sta, obs_tag):
        """get observed seismogram in array(ENZ,nt)"""
        if obs_tag not in g_sta:
            msg = f"{g_sta._v_name}: {obs_tag} does not exist, skip"
            raise KeyError(msg)
        obs = g_sta[obs_tag]

        # check data type
        obs_type = obs.attrs["type"]
        if obs_type != "VEL" and obs_type != "DISP":
            raise AssertionError(f"wrong types: obs({obs_type}), either VEL or DISP")

        # check channels' info
        obs_channels = obs.attrs["channels"]
        obs_Zchan = [cha for cha in obs_channels if cha["name"].decode()[-1] == "Z"]
        obs_Hchan = [cha for cha in obs_channels if cha["name"].decode()[-1] != "Z"]
        if len(obs_Zchan) != 1 or obs_Zchan[0]["dip"] != -90:
            msg = f"{g_sta._v_name} has problematic Z channel info: {obs_Zchan}, skip"
            raise AssertionError(msg)
        no_Hchan = False
        if len(obs_Hchan) == 0:
            msg = f"{g_sta._v_name} has no horizontal channels"
            warnings.warn(msg)
            no_Hchan = True
        else:
            assert len(obs_Hchan) == 2

        # rotate obs to ENZ
        # projection matrix: obs = proj * ENZ => ENZ = inv(proj) * obs
        data_nt = obs.attrs["npts"]
        obs_ENZ = np.zeros((3, data_nt), dtype=float)
        if no_Hchan:
            obs_ENZ[0:2, :] = 0.0
            obs_ENZ[2, :] = obs[0, :]  # only Z component
        else:
            assert obs.shape == obs_ENZ.shape
            proj_matrix = np.zeros((3, 3))
            for i in range(3):
                chan = obs_channels[i]
                sin_az = np.sin(np.deg2rad(chan["azimuth"]))
                cos_az = np.cos(np.deg2rad(chan["azimuth"]))
                sin_dip = np.sin(np.deg2rad(chan["dip"]))
                cos_dip = np.cos(np.deg2rad(chan["dip"]))
                # column vector = obs channel polarization
                proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
                proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
                proj_matrix[i, 2] = -sin_dip  # proj to Z (Up)
            # inverse projection matrix: ENZ = inv(proj) * obs
            inv_proj = np.linalg.inv(proj_matrix)
            obs_ENZ = np.dot(inv_proj, obs)

        return obs_ENZ, no_Hchan

    def __obsolete_get_syn_ENZ(self, g_sta, syn_tag):

        if syn_tag not in g_sta:
            msg = f"{g_sta._v_name}: {syn_tag} does not exist, skip"
            raise KeyError(msg)
        syn = g_sta[syn_tag]

        # check data type
        syn_type = syn.attrs["type"]
        if syn_type != "DISP":
            raise AssertionError(f"wrong types: syn({syn_type}) must be DISP")

        # check channels' info
        syn_channels = syn.attrs["channels"]
        assert len(syn_channels) == 3

        # rotate syn to ENZ
        data_nt = syn.attrs["npts"]
        data_fs = syn.attrs["sampling_rate"]
        data_dt = 1.0 / data_fs
        syn_ENZ = np.zeros((3, data_nt))
        syn_ENZ[:] = syn[:]
        proj_matrix = np.zeros((3, 3))
        for i in range(3):
            chan = syn_channels[i]
            sin_az = np.sin(np.deg2rad(chan["azimuth"]))
            cos_az = np.cos(np.deg2rad(chan["azimuth"]))
            sin_dip = np.sin(np.deg2rad(chan["dip"]))
            cos_dip = np.cos(np.deg2rad(chan["dip"]))
            # column vector = obs channel polarization
            proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
            proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
            proj_matrix[i, 2] = -sin_dip  # proj to Z
        # inverse projection matrix: ENZ = inv(proj) * obs
        inv_proj = np.linalg.inv(proj_matrix)
        syn_ENZ[:] = np.dot(inv_proj, syn_ENZ)

        # # apply time derivative (when data is vel)
        # if to_vel:
        #     npad = 100  # t_pad = 100 / fs, response drop to ~1% peak value
        #     nfft = scipy.fft.next_fast_len(data_nt + npad)
        #     freqs = np.fft.rfftfreq(nfft, d=data_dt)
        #     Fsyn *= 2j * np.pi * freqs
        #     syn_ENZ[:] = np.fft.irfft(Fsyn, nfft)[:, :data_nt]

        return syn_ENZ

    def __obsolete_extract_obs_syn_ENZ(self, g_sta, obs_tag, syn_tag, event_tau):
        if obs_tag not in g_sta:
            msg = f"{g_sta._v_name}: {obs_tag} does not exist, skip"
            raise KeyError(msg)
        if syn_tag not in g_sta:
            msg = f"{g_sta._v_name}: {syn_tag} does not exist, skip"
            raise KeyError(msg)

        obs = g_sta[obs_tag]
        syn = g_sta[syn_tag]

        # check data type
        obs_type = obs.attrs["type"]
        syn_type = syn.attrs["type"]
        case1 = obs_type == "VEL" and syn_type == "DISP"
        case2 = obs_type == "DISP" and syn_type == "DISP"
        if not (case1 or case2):
            raise AssertionError(f"wrong types: obs({obs_type}), syn({syn_type})")

        # check channels' info
        obs_channels = obs.attrs["channels"]
        syn_channels = syn.attrs["channels"]
        obs_Zchan = [cha for cha in obs_channels if cha["name"].decode()[-1] == "Z"]
        obs_Hchan = [cha for cha in obs_channels if cha["name"].decode()[-1] != "Z"]
        if len(obs_Zchan) != 1 or obs_Zchan[0]["dip"] != -90:
            msg = f"{g_sta._v_name} has problematic Z channel info: {obs_Zchan}, skip"
            raise AssertionError(msg)
        no_Hchan = False
        if len(obs_Hchan) == 0:
            msg = f"{g_sta._v_name} has no horizontal channels"
            warnings.warn(msg)
            no_Hchan = True
        else:
            assert len(obs_Hchan) == 2
        assert len(syn_channels) == 3

        # data_tb = obs.attrs["starttime"]
        data_fs = obs.attrs["sampling_rate"]
        data_dt = 1.0 / data_fs
        data_nt = obs.attrs["npts"]
        # data_te = data_tb + (data_nt - 1) * data_dt
        # data_times = np.arange(data_nt, dtype=float) * data_dt
        # noise_te = (first_arrtime - data_tb) - cfg_noise_before_first_arrival
        # noise_idx1 = int(noise_te * data_fs)  # 0:idx1 as noise

        # rotate obs to ENZ
        # projection matrix: obs = proj * ENZ => ENZ = inv(proj) * obs
        obs_ENZ = np.zeros((3, data_nt))
        if no_Hchan:
            obs_ENZ[2, :] = obs[0, :]  # only Z component
        else:
            assert obs.shape == obs_ENZ.shape
            proj_matrix = np.zeros((3, 3))
            for i in range(3):
                chan = obs_channels[i]
                sin_az = np.sin(np.deg2rad(chan["azimuth"]))
                cos_az = np.cos(np.deg2rad(chan["azimuth"]))
                sin_dip = np.sin(np.deg2rad(chan["dip"]))
                cos_dip = np.cos(np.deg2rad(chan["dip"]))
                # column vector = obs channel polarization
                proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
                proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
                proj_matrix[i, 2] = -sin_dip  # proj to Z (Up)
            # inverse projection matrix: ENZ = inv(proj) * obs
            inv_proj = np.linalg.inv(proj_matrix)
            obs_ENZ = np.dot(inv_proj, obs)

        # rotate syn to ENZ
        syn_ENZ = np.zeros((3, data_nt))
        syn_ENZ[:] = syn[:]
        proj_matrix = np.zeros((3, 3))
        for i in range(3):
            chan = syn_channels[i]
            sin_az = np.sin(np.deg2rad(chan["azimuth"]))
            cos_az = np.cos(np.deg2rad(chan["azimuth"]))
            sin_dip = np.sin(np.deg2rad(chan["dip"]))
            cos_dip = np.cos(np.deg2rad(chan["dip"]))
            # column vector = obs channel polarization
            proj_matrix[i, 0] = cos_dip * sin_az  # proj to E
            proj_matrix[i, 1] = cos_dip * cos_az  # proj to N
            proj_matrix[i, 2] = -sin_dip  # proj to Z
        # inverse projection matrix: ENZ = inv(proj) * obs
        inv_proj = np.linalg.inv(proj_matrix)
        syn_ENZ[:] = np.dot(inv_proj, syn_ENZ)

        # np.savez('extract', syn_ENZ=syn_ENZ)

        # apply source time function and/or time derivative
        if syn.attrs["is_grn"] or obs_type == "VEL":
            npad = int(5 * event_tau * data_fs)
            # npad = data_nt  #DEBUG
            nfft = scipy.fft.next_fast_len(data_nt + npad)
            freqs = np.fft.rfftfreq(nfft, d=data_dt)
            Fsyn = np.fft.rfft(syn_ENZ, nfft)
            if syn.attrs["is_grn"]:
                # source spectrum (moment-rate function)
                # print("event_tau=", event_tau)
                F_src = stf_gauss_spectrum(freqs, event_tau)
                Fsyn *= F_src
            if obs_type == "VEL":
                Fsyn *= 2j * np.pi * freqs
            syn_ENZ[:] = np.fft.irfft(Fsyn, nfft)[:, :data_nt]

        return obs_ENZ, syn_ENZ, no_Hchan

    def __obsolete_measure_adj(self):
        """
        calculate adjoint sources (dchi_du)

        Parameters
        ----------
        misfit_type:
          cc0: zero-lag cross-correlation misfit
          ccdt: cross-correlation traveltime
          phcc: phase correlation (for ambient noise correlogram)
        weight_param : SNR, CCmax, CC0, cc_tshift

        Notes
        -----
        chi : misfit value (normalized zero-lag correlation coef.)
        u : synthetic waveform

        """
        h5f = pt.open_file(self.h5_path, "r+")

        config = h5f.root._v_attrs["config"]
        obs_tag = config["data"]["tag"]
        syn_tag = config["syn"]["tag"]
        cfg_noise_before_first_arrival = config["data"][
            "time_before_first_arrival_as_noise"
        ]
        cci_dt = config["misfit"]["cc_deltat"]
        misfit_type = config["misfit"]["type"]
        weight_param = config["misfit"]["window_weight"]
        adj_tag = config["adj"]["tag"]
        adj_band_code = config["adj"]["band_code"]

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])
        event_tau = event["tau"]

        if "/channel" not in h5f:
            msg = '"/channel" not existing, run read_channel_file first!'
            raise KeyError(msg)
        tbl_chan = h5f.get_node("/channel")

        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        if "/window" not in h5f:
            msg = '"/waveform" not existing, run setup_windows first!'
            raise KeyError(msg)
        tbl_win = h5f.get_node("/window")

        # for storing adjoint source
        pt_atom = pt.Atom.from_dtype(np.dtype(np.float32))
        pt_filters = pt.Filters(complevel=3, complib="zlib")

        # loop each station
        for g_sta in g_wav:
            # print(f"measure_adj for {g_sta._v_name}")

            attrs = g_sta._v_attrs
            net = attrs["network"]
            sta = attrs["station"]
            stnm = f"{net}_{sta}"
            # loc = attrs['location']
            first_arrtime = attrs["first_arrtime"]

            if obs_tag not in g_sta:
                msg = f"{obs_tag} does not exist under {g_sta._v_file.filename}:{g_sta._v_pathname}, skip"
                warnings.warn(msg)
                continue
            if syn_tag not in g_sta:
                msg = f"{syn_tag} does not exist under {g_sta._v_file.filename}:{g_sta._v_pathname}, skip"
                warnings.warn(msg)
                continue

            obs = g_sta[obs_tag]
            syn = g_sta[syn_tag]

            sta_attrs = {k: g_sta._v_attrs[k] for k in g_sta._v_attrs._f_list("user")}
            obs_attrs = {k: obs.attrs[k] for k in obs.attrs._f_list("user")}
            syn_attrs = {k: syn.attrs[k] for k in syn.attrs._f_list("user")}
            try:
                # obs_ENZ, syn_ENZ, no_Hchan = self._extract_obs_syn_ENZ(
                #     g_sta, obs_tag, syn_tag, event["tau"]
                # )
                obs_ENZ, syn_ENZ, no_Hchan = _extract_obs_syn_ENZ(
                    obs[:],
                    obs_attrs,
                    syn[:],
                    syn_attrs,
                    # sta_attrs, obs[:], obs_attrs, syn[:], syn_attrs
                )
                # debug
                # np.savez(f'{net}_{sta}_old.npz', obs_ENZ=obs_ENZ, syn_ENZ=syn_ENZ, no_Hchan=no_Hchan)
            except Exception as e:
                msg = f"{stnm}: failed to get obs_ENZ,syn_ENZ ({e})"
                warnings.warn(msg)
                continue

            obs_type = obs_attrs["type"]
            syn_type = syn_attrs["type"]
            if obs_type == "VEL" and syn_type == "DISP":
                syn_to_vel = True
            elif obs_type == "DISP" and syn_type == "DISP":
                syn_to_vel = False
            else:
                msg = f"{stnm}: unknown obs/syn types ({obs_type}/{syn_type})"
                raise ValueError(msg)
            conv_stf = False
            if syn_attrs["is_grn"]:
                conv_stf = True
                assert syn_attrs["origin_time"] == event_t0

            data_tb = obs.attrs["starttime"]
            data_fs = obs.attrs["sampling_rate"]
            data_dt = 1.0 / data_fs
            data_nt = obs.attrs["npts"]
            data_te = data_tb + (data_nt - 1) * data_dt
            data_times = np.arange(data_nt, dtype=float) * data_dt
            noise_te = (first_arrtime - data_tb) - cfg_noise_before_first_arrival
            # noise_idx1 = int(noise_te * data_fs)  # 0:idx1 as noise

            # pre-filter parameters
            pre_filters = obs.attrs["filter"]
            # data_butter_N = obs.attrs["filter"]["N"]
            # data_butter_Wn = obs.attrs["filter"]["Wn"]

            # sum of adjoint sources from all windows
            dchi_du = np.zeros((3, data_nt))

            # loop each window
            for win in tbl_win.where(f'(network == b"{net}") & (station == b"{sta}")'):
                # print(
                #     "window: ",
                #     win["id"],
                #     win["type"],
                #     win["phase"],
                #     win["cmpnm"],
                #     win["butter_Wn"],
                # )
                # component
                cmpnm = win["cmpnm"].decode()
                if no_Hchan and cmpnm != "Z":
                    warnings.warn(
                        "skip non-vertical window since only vertical data present"
                    )
                    win["valid"] = False
                    win.update()
                    continue
                cmpaz = win["cmpaz"]
                cmpdip = win["cmpdip"]

                # filter
                butter_N = win["butter_N"]
                butter_Wn = win["butter_Wn"]
                # fft
                npad = int(2 * data_fs / min(butter_Wn))
                nfft = scipy.fft.next_fast_len(data_nt + npad)
                freqs = np.fft.rfftfreq(nfft, d=data_dt)
                # window filter
                sos = scipy.signal.butter(
                    butter_N, butter_Wn, "bandpass", fs=data_fs, output="sos"
                )
                _, filter_h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
                Fw = abs(filter_h)  # zero-phase filter response
                # pre-filter applied to obs & syn
                filter_h = np.ones_like(freqs)
                for filter in pre_filters:
                    sos = scipy.signal.butter(
                        filter["N"],
                        filter["Wn"],
                        filter["type"],
                        fs=data_fs,
                        output="sos",
                    )
                    _, h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
                    filter_h *= abs(h)
                Fd = abs(filter_h)
                # time derivative in frequency domain
                Ft = 2j * np.pi * freqs
                # time window
                win_tb = UTCDateTime(win["starttime"])
                win_te = UTCDateTime(win["endtime"])
                win_taper = win["taper"]

                # window function
                b = win_tb - data_tb
                e = win_te - data_tb
                win_len = e - b
                width = win_len * min(win_taper, 0.4)
                win_c = [b, b + width, e - width, e]
                win_func = cosine_sac_taper(data_times, win_c)

                # component
                if cmpnm in ["Z", "R", "T"]:
                    sin_az = np.sin(np.deg2rad(cmpaz))
                    cos_az = np.cos(np.deg2rad(cmpaz))
                    sin_dip = np.sin(np.deg2rad(cmpdip))
                    cos_dip = np.cos(np.deg2rad(cmpdip))
                    n = np.array(
                        [
                            [cos_dip * sin_az],  # cos(E, comp)
                            [cos_dip * cos_az],  # cos(N, comp)
                            [-sin_dip],  # cos(Z, comp)
                        ]
                    )
                    proj_matrix = np.dot(n, n.transpose())
                elif cmpnm == "H":  # horizontal vector 2d
                    proj_matrix = np.identity(3)
                    proj_matrix[2, 2] = 0.0  # zero Z component
                elif cmpnm == "F":  # full 3d vector
                    proj_matrix = np.identity(3)
                else:
                    msg = f"unrecognized component code ({cmpnm}), SKIP"
                    warnings.warn(msg)
                    continue

                # DEBUG
                # np.save("obs_ENZ.npy", obs_ENZ)
                # np.save("syn_ENZ.npy", syn_ENZ)
                # np.save("Fd.npy", Fd)
                # np.save("Fw.npy", Fw)
                # np.save("win_func.npy", win_func)
                # np.save("proj_matrix.npy", proj_matrix)
                # print(data_nt, nfft, data_fs)

                # Fd: pre-filter used during preprocessing (e.g. rmresp)
                # Fw: bandpass filter used on this window
                # w: window function
                #
                # obs = w*(Fw*Fd*d), syn = w*(Fw*Fd*u)

                # apply filter to obs, syn
                # obs = Fw * (Fd * d), obs_ENZ = Fd * d, d = disp. or vel.
                obs_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ, nfft), nfft)[
                    :, :data_nt
                ]
                # syn = Fw * (Fd * [Ft] * [S] * u), (syn_ENZ = Fd * u)
                f_syn = np.fft.rfft(syn_ENZ, nfft)
                if syn_to_vel:
                    f_syn *= Ft
                if conv_stf:
                    f_src = stf_gauss_spectrum(freqs, event_tau)
                    f_syn *= f_src
                syn_filt = np.fft.irfft(Fw * f_syn, nfft)[:, :data_nt]
                # # syn = Fw * (Fd * [Ft] * u), (syn_ENZ = Fd * [Ft] *[S] * u)
                # syn_filt = np.fft.irfft(Fw * np.fft.rfft(syn_ENZ, nfft), nfft)[
                #     :, :data_nt
                # ]
                # noise
                taper_width = 0.5 / min(butter_Wn)
                noise_win = cosine_sac_taper(
                    data_times, [0, taper_width, noise_te - taper_width, noise_te]
                )
                noise_idx = (data_times > taper_width) & (
                    data_times < noise_te - taper_width
                )
                # # cut noise window from obs_ENZ before applying bandpass filter Fw
                # noise_filt = np.fft.irfft(
                #     Fw * np.fft.rfft(obs_ENZ * noise_win, nfft), nfft
                # )[:, :data_nt]

                # # DEBUG
                # for i in range(3):
                #     plt.subplot(311 + i)
                #     plt.plot(
                #         data_times, obs[i, :], "k", data_times, obs_filt[i, :], "r", data_times, noise_filt[i,:], 'b'
                #     )
                # plt.show()
                # for i in range(3):
                #     plt.subplot(311 + i)
                #     plt.plot(
                #         data_times, syn[i, :], "k", data_times, syn_filt[i, :], "r"
                #     )
                # plt.show()

                # apply window and projection
                # obs: w*(Fw*Fd*d)
                obs_filt_win = np.dot(proj_matrix, obs_filt) * win_func
                # syn: w*(Fw*Fd*u)
                syn_filt_win = np.dot(proj_matrix, syn_filt) * win_func
                # noise
                noise_filt_win = np.dot(proj_matrix, obs_filt) * noise_win

                # np.savez('measure_adj', obs_ENZ=obs_ENZ, syn_ENZ=syn_ENZ, syn_filt=syn_filt, win_func=win_func, obs_filt_win=obs_filt_win, syn_filt_win=syn_filt_win)

                # for i in range(3):
                #     plt.subplot(311 + i)
                #     plt.plot(
                #         data_times, syn[i, :], "k", data_times, syn_filt[i, :], "r"
                #     )
                # plt.show()

                # DEBUG
                # diff = obs_ENZ_win - syn_ENZ_win
                # for i in range(3):
                #  plt.subplot(311+i)
                #  plt.plot(syn_times, obs_ENZ_win[i,:], 'k')
                #  plt.plot(syn_times, syn_ENZ_win[i,:], 'r')
                #  plt.plot(syn_times, diff[i,:], 'c')
                # plt.show()

                # ------ measure SNR (based on maximum amplitude)
                Amax_obs = np.sqrt(np.max(np.sum(obs_filt_win**2, axis=0)))
                Amax_syn = np.sqrt(np.max(np.sum(syn_filt_win**2, axis=0)))
                Amax_noise = np.sqrt(
                    np.max(np.sum(noise_filt_win[:, noise_idx] ** 2, axis=0))
                )
                if Amax_obs == 0:  # bad record
                    msg = f"empty obs trace ({win}), skip"
                    warnings.warn(msg)
                    continue
                if Amax_noise == 0:
                    # could occure when the data begin time is too close to the first arrival
                    msg = f"empty noise trace ({win}), SKIP."
                    warnings.warn(msg)
                    continue
                snr = 20.0 * np.log10(Amax_obs / Amax_noise)

                # ------ measure CC time shift (between w*F*d and w*F*u)
                obs_norm = np.sqrt(np.sum(obs_filt_win**2))
                syn_norm = np.sqrt(np.sum(syn_filt_win**2))
                # window normalization factor (without dt)
                Nw = obs_norm * syn_norm
                # NOTE the order (obs,syn) is important. The positive time on
                # CC means shifting syn in the positive time direction to match
                # the observed obs, and vice verser.
                # [-(nt-1), nt) * dt
                cc = np.zeros(2 * data_nt - 1)
                for i in range(3):
                    cc += signal.fftconvolve(
                        obs_filt_win[i, :], syn_filt_win[i, ::-1], "full"
                    )
                cc /= Nw
                # -- zero-lag cc coeff.
                CC0 = cc[data_nt - 1]  # the n-th point corresponds to zero lag time
                AR0 = CC0 * syn_norm / obs_norm  # amplitude ratio syn/obs
                # DEBUG
                # print(CC0, np.sum(obs_filt_win * syn_filt_win)/obs_norm/syn_norm)
                # -- interpolate cc to finer time samples
                CC_shift_range = win_len / 2.0  # TODO: more reasonable choice?
                if data_dt < cci_dt:
                    msg = f"data_dt({data_dt}) < cci_dt({cci_dt})"
                    warnings.warn(msg)
                    cci_dt = data_dt
                ncci = int(CC_shift_range / cci_dt)
                cci_t0 = -ncci * cci_dt
                cc_t0 = -(data_nt - 1) * data_dt  # begin time in cc
                cci = lanczos_interpolation(
                    cc, cc_t0, data_dt, cci_t0, cci_dt, 2 * ncci + 1, a=20
                )
                # time shift at the maximum correlation
                imax = np.argmax(cci)
                CC_time_shift = cci_t0 + imax * cci_dt
                CCmax = cci[imax]
                ARmax = CCmax * syn_norm / obs_norm  # amplitude ratio: syn/obs

                # ------ window weighting based on SNR and misfit
                weight = win["weight"]
                if "SNR" in weight_param:
                    weight *= cosine_taper(snr, weight_param["SNR"])
                if "CCmax" in weight_param:
                    weight *= cosine_taper(CCmax, weight_param["CCmax"])
                if "CC0" in weight_param:
                    weight *= cosine_taper(CC0, weight_param["CC0"])
                if "cc_tshift" in weight_param:
                    weight *= cosine_taper(CC_time_shift, weight_param["cc_tshift"])

                # ------ measure adjoint source
                Aw = CC0 * obs_norm / syn_norm  # window amplitude raito
                if misfit_type == "cc0":
                    #
                    # Goodness-of-fit: (zero-lag normalized cross-correlation coefficient)
                    #   cc(syn, obs) = sum(syn * obs) / sqrt(sum(syn**2) * sum(obs**2))
                    #
                    # obs = w * (Fw * Fd * P * d),        d = displacement or velocity.
                    # syn = w * (Fw * Fd * P * u),        if d = disp.
                    #       w * (Fw * Fd * Ft * P * u),   if d = vel.
                    #
                    #  u[3,nt]: synthetic seismogram of displacement
                    #     , replaced by S * u, when u is Green's function (S: source wavelet)
                    #  d[3,nt]: recorded ground displacement or velocity
                    #  P[3,3]: the projection matrix (real valued), e.g. project to radial/tangential component
                    #  Fd[nfft]: bandpass filter applied during pre-processing (e.g. removing instrument response, down-sampling)
                    #  Fw[nfft]: bandpass filter applied for the measurment time window
                    #  Ft[nfft]: time derivate (2j * pi * freqs in frequency domian)
                    #  note: F * x = irfft(F * rfft(x, nfft), nfft)[:, :nt]
                    #  w[nt]: time window function
                    #
                    # adjoint source:
                    #   Dcc/Du[3,nt] = conj(Fw * Fd * [Ft]) * transpose(P) * w * (obs - Aw * syn) / obs_norm / syn_norm,
                    #
                    #   where,
                    #           [Ft] is only applied when data is velocity.
                    #           obs_norm = sqrt(sum(obs**2))
                    #           syn_norm = sqrt(sum(syn**2))
                    #           Aw = sum(syn * obs) / sum(syn**2)
                    #
                    #   cc(u+du) - cc(u) = sum(Dcc/Du * du) + O(du**2)
                    #
                    dchiw_du = np.dot(
                        np.transpose(proj_matrix),
                        win_func * (obs_filt_win - Aw * syn_filt_win) / Nw,
                    )
                    conj_F = np.conjugate(Fw * Fd, dtype="complex")
                    if obs.attrs["type"] == "VEL":
                        conj_F *= np.conjugate(Ft)
                    dchiw_du = np.fft.irfft(conj_F * np.fft.rfft(dchiw_du, nfft), nfft)[
                        :, :data_nt
                    ]
                    # np.save("dchiw_du.npy", dchiw_du)

                # elif misfit_type == "phcc":
                #     # misfit: weighted phase correlation
                #     norm_N = np.sum(obs_filt_win**2)
                #     fft_wFu = np.fft.rfft(syn_filt_win)
                #     fft_wFd = np.fft.rfft(obs_filt_win)
                #     abs_fft_wFu = (
                #         np.sum(np.abs(fft_wFu) ** 2, axis=0, keepdims=True) ** 0.5
                #     )
                #     abs_fft_wFd = (
                #         np.sum(np.abs(fft_wFd) ** 2, axis=0, keepdims=True) ** 0.5
                #     )
                #     abs_fft_wFu_wl = np.copy(abs_fft_wFu)
                #     wl_thred = 0.01 * np.max(
                #         abs_fft_wFu
                #     )  # FIXME water-level coef. 0.01 is hard coded here!
                #     abs_fft_wFu_wl[abs_fft_wFu < wl_thred] = wl_thred
                #     trans_H = abs_fft_wFd / abs_fft_wFu_wl
                #     dchiw_du = (
                #         trans_H
                #         / norm_N
                #         * (
                #             fft_wFd
                #             - fft_wFu
                #             * np.real(
                #                 np.sum(
                #                     fft_wFd * np.conj(fft_wFu), axis=0, keepdims=True
                #                 )
                #             )
                #             / abs_fft_wFu_wl**2
                #         )
                #     )
                #     dchiw_du = np.fft.irfft(dchiw_du) * win_func
                #     dchiw_du = signal.filtfilt(filter_b, filter_a, dchiw_du[:, ::-1])
                #     dchiw_du = dchiw_du[:, ::-1]
                #     # phase correlation value
                #     HwFu = np.fft.irfft(trans_H * fft_wFu)
                #     phcc = np.sum(obs_filt_win * HwFu) / N
                # elif misfit_type == "ccdt":
                #     # misfit: cross-correlation traveltime difference -1*|ccdtau|^2
                #     # adjoint source: dchiw_du = ccdt * dtau/du
                #     # dchiw_du = conj(F)*w*d/dt(w*F*u)/N,
                #     # , where N = norm(d/dt(w*F*u))
                #     # -- dchiw_du
                #     # NOTE: *dt is put back to Nw
                #     dwFu_dt = np.fft.irfft(
                #         np.fft.rfft(syn_filt_win) * 2.0j * np.pi * syn_freq, n=syn_nt
                #     )
                #     w_dwFu_dt = win_func * dwFu_dt
                #     # apply conj(F), for two-pass filter (zero phase) conj(F) = F
                #     conjF_w_dwFu_dt = signal.filtfilt(
                #         filter_b, filter_a, w_dwFu_dt[:, ::-1]
                #     )
                #     conjF_w_dwFu_dt = conjF_w_dwFu_dt[:, ::-1]
                #     #
                #     dchiw_du = -1 * CC_time_shift * conjF_w_dwFu_dt / np.sum(dwFu_dt**2)
                #     # DEBUG
                #     # print(CC_time_shift)
                #     # plt.subplot(411)
                #     # plt.plot(syn_times, syn_filt_win[1,:])
                #     # plt.subplot(412)
                #     # plt.plot(syn_times, dwFu_dt[1,:])
                #     # plt.subplot(413)
                #     # plt.plot(syn_times, conjF_w_dwFu_dt[1,:])
                #     # plt.subplot(414)
                #     # plt.plot(syn_times, dchiw_du[1,:])
                #     # plt.show()
                else:
                    msg = f"{g_sta._v_name}: unknown misfit type ({misfit_type})"
                    raise ValueError(msg)

                # DEBUG
                # for i in range(3):
                #  plt.subplot(311+i)
                #  plt.plot(syn_times, dchiw_du1[i,:], 'k')
                #  plt.plot(syn_times, dchiw_du[i,:], 'r')
                # plt.show()

                # add into total dchi_du
                dchi_du += weight * dchiw_du

                # -- dchiw_dg = conj(S) * dchiw_du
                # dchiw_dg = np.fft.irfft(np.conjugate(F_src) * np.fft.rfft(dchiw_du), syn_nt)
                # add into total dchi_dg
                # dchi_dg += weight * dchiw_dg
                # DEBUG
                # for i in range(3):
                #  plt.subplot(311+i)
                #  plt.plot(syn_times, dchiw_du[i,:], 'k')
                #  plt.plot(syn_times, dchiw_dg[i,:], 'r')
                # plt.show()

                # DEBUG
                if _DEBUG:
                    print(f"cc_time_shift = {CC_time_shift}")
                    print(f"cc_max = {CCmax}")
                    print(f"SNR = {snr}")
                    print(f"CC0 = {CC0}")
                    print(f"Aw = {Aw}")
                    phase_shift = np.exp(-2j * np.pi * freqs * CC_time_shift)
                    Fsyn = np.fft.rfft(syn_filt_win, nfft) * phase_shift
                    syn_filt_win_shift = np.fft.irfft(Fsyn, nfft)[:, :data_nt]
                    scale_dcc = np.max(np.abs(syn_filt_win)) / np.max(np.abs(dchiw_du))
                    print(f"{scale_dcc=}")
                    du = obs_filt_win - Aw * syn_filt_win
                    for i in range(3):
                        plt.subplot(311 + i)
                        plt.plot(
                            data_times,
                            obs_filt_win[i, :],
                            "k",
                            # data_times, obs_ENZ[i, :], "k--",
                            data_times,
                            syn_filt_win[i, :],
                            "r",
                            # data_times, syn_ENZ[i, :], "r--",
                            # data_times, syn_filt_win_shift[i, :], "c",
                            # data_times[:noise_idx1], noise_filt_win[i,:], "b",
                            # data_times, scale_dcc * dchiw_du[i,:], "y"
                            data_times,
                            1.0e6 * du[i, :],
                            "y",
                        )
                        if i == 0:
                            plt.title(f"{win['id']}")
                    plt.show()

                win["cc0"] = CC0
                win["cc_time_shift"] = CC_time_shift
                win["cc_max"] = CCmax
                win["amp_ratio_cc0"] = AR0
                win["amp_ratio_ccmax"] = ARmax
                win["noise_maxamp"] = Amax_noise
                win["obs_maxamp"] = Amax_obs
                win["syn_maxamp"] = Amax_syn
                win["SNR"] = snr
                win["weight"] = weight
                win["valid"] = True
                win["proj_matrix"] = proj_matrix
                win.update()

            # DEBUG
            # wins = tbl_win.read_where(f'(network == b"{net}") & (station == b"{sta}")')
            # np.savez(f'{net}_{sta}_wins.npz', wins=wins)

            # store adjoint source for this station, e.g. /NET_STA/ADJ_DISP[0:nchan, 0:npts]
            if adj_tag in g_sta:
                msg = f"{adj_tag} exists in {g_sta._v_name}, overwrite!"
                warnings.warn(msg)
                h5f.remove_node(g_sta, adj_tag)
            shape = (3, data_nt)
            ca = h5f.create_carray(g_sta, adj_tag, pt_atom, shape, filters=pt_filters)
            ca[:] = dchi_du
            dtype = [("name", "S3"), ("azimuth", float), ("dip", float)]
            channels = [
                (f"{adj_band_code}E", 90, 0),
                (f"{adj_band_code}N", 0, 0),
                (f"{adj_band_code}Z", 0, -90),
            ]
            ca.attrs["channels"] = np.array(channels, dtype=dtype)
            ca.attrs["starttime"] = data_tb
            ca.attrs["sampling_rate"] = data_fs
            ca.attrs["npts"] = data_nt
            ca.attrs["solver_tb"] = syn.attrs["solver_tb"]
            ca.attrs["solver_dt"] = syn.attrs["solver_dt"]
            ca.attrs["solver_nt"] = syn.attrs["solver_nt"]
            ca.attrs["type"] = syn.attrs["type"]
            ca.attrs["misfit"] = misfit_type

            # DEBUG
            # for i in range(3):
            #  plt.subplot(311+i)
            #  plt.plot(syn_times, dchi_du[i,:], 'k')
            #  plt.plot(syn_times, dchi_dg[i,:], 'r')
            # plt.show()

        h5f.close()

    class FileAccess(multiprocessing.Process):

        def __init__(self, h5_path, read_queue, result_queues, write_queue, shutdown):
            self.h5_path = h5_path
            self.read_queue = read_queue
            self.result_queues = result_queues
            self.write_queue = write_queue
            self.shutdown = shutdown
            self.block_period = 0.01
            super().__init__()

        def run(self):
            self.h5_file = pt.open_file(self.h5_path, "r+")
            stop_loop = False
            end_of_read_queue = False
            end_of_write_queue = False
            while True:
                if self.shutdown.poll():
                    stop_loop = True
                if stop_loop and end_of_read_queue and end_of_write_queue:
                    break

                try:
                    item = self.read_queue.get(True, self.block_period)
                    if item is None:
                        end_of_read_queue = True
                    else:
                        proc_num, net, sta, obs_tag, syn_tag = item
                        data = self.read_data(net, sta, obs_tag, syn_tag)
                        result_queue = self.result_queues[proc_num]
                        result_queue.put(data)
                except queue.Empty:
                    pass

                try:
                    item = self.write_queue.get(True, self.block_period)
                    if item is None:
                        end_of_write_queue = True
                    else:
                        net, sta, adj_tag, adj_src, adj_attrs, inds, wins = item
                        self.write_data(
                            net, sta, adj_tag, adj_src, adj_attrs, inds, wins
                        )
                except queue.Empty:
                    pass
            # close the HDF5 file before shutting down
            self.h5_file.close()

        def read_data(self, net, sta, obs_tag, syn_tag):
            try:
                gsta = self.h5_file.get_node(f"/waveform/{net}_{sta}")
                tbl_win = self.h5_file.get_node("/window")

                sta_attrs = {k: gsta._v_attrs[k] for k in gsta._v_attrs._f_list("user")}

                obs = gsta[obs_tag]
                obs_data = obs[:]
                obs_attrs = {k: obs.attrs[k] for k in obs.attrs._f_list("user")}

                syn = gsta[syn_tag]
                syn_data = syn[:]
                syn_attrs = {k: syn.attrs[k] for k in syn.attrs._f_list("user")}

                inds = tbl_win.get_where_list(
                    f'(network == b"{net}") & (station == b"{sta}")'
                )
                wins = tbl_win[inds]
            except Exception as e:
                msg = f"{net}_{sta}: failed to get syn/obs data"
                warnings.warn(msg)
                return None

            return (
                sta_attrs,
                obs_data,
                obs_attrs,
                syn_data,
                syn_attrs,
                inds,
                wins,
            )

        def write_data(self, net, sta, adj_tag, adj_src, adj_attrs, inds, wins):
            gsta = self.h5_file.get_node(f"/waveform/{net}_{sta}")
            tbl_win = self.h5_file.get_node("/window")

            pt_atom = pt.Atom.from_dtype(np.dtype(np.float32))
            pt_filters = pt.Filters(complevel=3, complib="zlib")

            if adj_tag in gsta:
                msg = f"{net}_{sta}: {adj_tag} exists, overwrite!"
                warnings.warn(msg)
                self.h5_file.remove_node(gsta, adj_tag)

            ca = self.h5_file.create_carray(
                gsta, adj_tag, pt_atom, adj_src.shape, filters=pt_filters
            )
            ca[:] = adj_src
            for k in adj_attrs:
                ca.attrs[k] = adj_attrs[k]

            tbl_win[inds] = wins
            tbl_win.flush()

    class DataProcessor(multiprocessing.Process):

        def __init__(
            self,
            read_queue,
            result_queue,
            write_queue,
            proc_num,
            net_sta_list,
            config,
            event_t0,
            event_tau,
        ):
            self.read_queue = read_queue
            self.result_queue = result_queue
            self.write_queue = write_queue
            self.proc_num = proc_num
            self.net_sta_list = net_sta_list
            self.config = config
            self.event_t0 = event_t0
            self.event_tau = event_tau
            super().__init__()

        def run(self):
            obs_tag = self.config["data"]["tag"]
            syn_tag = self.config["syn"]["tag"]
            cfg_noise_before_first_arrival = self.config["data"][
                "time_before_first_arrival_as_noise"
            ]
            cci_dt = self.config["misfit"]["cc_deltat"]
            misfit_type = self.config["misfit"]["type"]
            weight_param = self.config["misfit"]["window_weight"]
            adj_tag = self.config["adj"]["tag"]
            adj_band_code = self.config["adj"]["band_code"]
            event_t0 = self.event_t0
            event_tau = self.event_tau

            for net, sta in self.net_sta_list:
                self.read_queue.put((self.proc_num, net, sta, obs_tag, syn_tag))
                ret = self.result_queue.get()
                if not ret:
                    continue
                sta_attrs, obs_data, obs_attrs, syn_data, syn_attrs, inds, wins = ret

                ret = _measure_adj_one_sta(
                    sta_attrs,
                    obs_data,
                    obs_attrs,
                    syn_data,
                    syn_attrs,
                    wins,
                    event_t0,
                    event_tau,
                    cfg_noise_before_first_arrival,
                    cci_dt,
                    misfit_type,
                    weight_param,
                    adj_band_code,
                )
                if ret:
                    adj_src, adj_attrs, wins = ret
                    self.write_queue.put(
                        (net, sta, adj_tag, adj_src, adj_attrs, inds, wins)
                    )

    def measure_adj_multiprocess(self, nproc=5):
        """
        calculate adjoint sources (dchi_du)
        """
        assert nproc >= 1

        h5f = pt.open_file(self.h5_path, "r+")

        config = h5f.root._v_attrs["config"]

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])
        event_tau = event["tau"]

        if "/channel" not in h5f:
            msg = '"/channel" not existing, run read_channel_file first!'
            raise KeyError(msg)
        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")
        if "/window" not in h5f:
            msg = '"/waveform" not existing, run setup_windows first!'
            raise KeyError(msg)

        net_sta_list = [
            (gsta._v_attrs["network"], gsta._v_attrs["station"]) for gsta in g_wav
        ]

        h5f.close()

        multiprocessing.set_start_method("spawn")

        read_queue = multiprocessing.Queue()
        write_queue = multiprocessing.Queue()
        shutdown_recv, shutdown_send = multiprocessing.Pipe(False)
        result_queues = [multiprocessing.Queue() for i in range(nproc)]
        file_access = self.FileAccess(
            self.h5_path, read_queue, result_queues, write_queue, shutdown_recv
        )
        file_access.start()
        processors = []
        for i in range(nproc):
            result_queue = result_queues[i]
            processor = self.DataProcessor(
                read_queue,
                result_queue,
                write_queue,
                i,
                net_sta_list[i::nproc],
                config,
                event_t0,
                event_tau,
            )
            processors.append(processor)
        for processor in processors:
            processor.start()
        for processor in processors:
            processor.join()
        # shut down the FileAccess instance
        read_queue.put(None)
        write_queue.put(None)
        shutdown_send.send(0)
        file_access.join()

    def plot_seismogram_1comp(
        self,
        # savefig=False,
        # out_dir="plot",
        win_id,  # e.g. p,P_Z_30-100sec
        azbin=20,
        max_ntrace_per_bin=30,  # maximum traces to plot for each azimuthal bin
        begin_time=0,
        end_time=0,
        clip_ratio=1.5,
        min_CC0=None,
        min_CCmax=None,
        min_SNR=None,
        dist_lim=None,
        # plot_az0=0,
        # plot_adj=False,  # whether plot adjoint source
        align_time=False,  # whether align the phase according to cc time shift
        out_fig="seis.pdf",
        title_prefix=None,
    ):
        """
        Plot record section in azimuthal bins

        Parameters
        ----------
        azbin: azimuthal bin size

        begin/end_time: time range that is added to the automatically determined
          plot time range. See below.

        clip: do not plot waveform with amplitudes larger than
          <clip>*max_amplitude_in_select_time_window

        Notes
        -----
          The time for plot is reduced time relative to origin time + rayp*dist

          linear regression is done to find the average rayp of misfit windows
          and the plot begin/end time is found that can include all misfit windows
          in the plot.

        """
        # ------ check parameters
        # plot_time = np.array([begin_time, end_time])

        # in case reverse the distance axis
        # plot_flip = -1
        plot_flip = 1

        plot_azbin = float(azbin)
        if plot_azbin <= 0:
            raise Exception("plot_azbin(%f) should be larger than 0.0" % (plot_azbin))
        try:
            plot_max_ntrace_per_bin = int(max_ntrace_per_bin)
            if plot_max_ntrace_per_bin < 1:
                plot_max_ntrace_per_bin = None
        except:
            plot_max_ntrace_per_bin = None
        plot_window_id = win_id
        plot_SNR = np.array(min_SNR)
        plot_CC0 = np.array(min_CC0)
        plot_CCmax = np.array(min_CCmax)
        plot_distlim = np.array(dist_lim)

        plot_clip = float(clip_ratio)
        if plot_clip < 1.0:
            raise Exception("clip_ratio(%f) should be larger than 1.0" % (plot_clip))

        h5f = pt.open_file(self.h5_path, "r")

        # config
        config = h5f.root._v_attrs["config"]
        taup_model = TauPyModel(model=config["taup_model"])
        obs_tag = config["data"]["tag"]
        syn_tag = config["syn"]["tag"]

        # event info
        if "/source" not in h5f:
            msg = '"/source" does not exist, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]

        if "/channel" not in h5f:
            msg = '"/channel" does not exist, run read_channel_file first!'
            raise KeyError(msg)
        tbl_chan = h5f.get_node("/channel")

        if "/waveform" not in h5f:
            msg = '"/waveform" does not exist, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        if "/window" not in h5f:
            msg = '"/waveform" does not exist, run setup_windows first!'
            raise KeyError(msg)
        tbl_win = h5f.get_node("/window")

        # event info
        evt0 = UTCDateTime(event["t0"])
        evtau = event["tau"]
        evla = event["latitude"]
        evlo = event["longitude"]
        evdp = event["depth"]
        evnm = event["id"].decode()
        # evdp has to be >=0 otherwise taup would crash
        if evdp < 0.0:
            evdp = 0.0
        mt = event["mt_rtp"]
        Mrr = mt[0][0]
        Mtt = mt[1][1]
        Mpp = mt[2][2]
        Mrt = mt[0][1]
        Mrp = mt[0][2]
        Mtp = mt[1][2]
        focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]

        # get windows
        cmpnm_select = ["Z", "R", "T"]
        windows = [
            (win, g_wav[stnm])
            for win in tbl_win.read_where(f'id == b"{win_id}"')
            if (stnm := f'{win["network"].decode()}_{win["station"].decode()}') in g_wav
            and win["cmpnm"].decode() in cmpnm_select
        ]
        if not windows:
            warnings.warn("No data to plot!")
            return
        stla_all = np.array([sta._v_attrs["latitude"] for _, sta in windows])
        stlo_all = np.array([sta._v_attrs["longitude"] for _, sta in windows])
        dist_all = np.array([sta._v_attrs["dist_degree"] for _, sta in windows])
        weight_all = np.array([win["weight"] for win, _ in windows])
        ccdt_all = np.array([win["cc_time_shift"] for win, _ in windows])
        cc0_all = np.array([win["cc0"] for win, _ in windows])
        winb_all = np.array(
            [UTCDateTime(win["starttime"]) - evt0 for win, _ in windows]
        )
        wine_all = np.array([UTCDateTime(win["endtime"]) - evt0 for win, _ in windows])
        winc_all = (winb_all + wine_all) / 2.0
        # get average moveout of the window center
        # linear regression tc = dist*rayp + tb
        A = np.vstack([dist_all, np.ones(len(dist_all))]).T
        plot_rayp, plot_c = np.linalg.lstsq(A, winc_all, rcond=None)[0]
        # round to integer
        plot_rayp = np.round(plot_rayp)

        # # map configuration
        # parallels = np.arange(-90.0, 90, 10.0)
        # meridians = np.arange(0.0, 360, 10.0)
        # map_lat0, map_lon0, map_width, map_height = centerMap(
        #     [*stla_all, evla], [*stlo_all, evlo], 1.1
        # )

        # ------ calculate traveltime curves (only for body wave)
        phases = set([w["phase"].decode() for w, _ in windows if w["type"] == b"body"])
        phase_list = [a for p in phases for a in p.split(",")]
        if windows[0][0]["type"] == b"Pwin":
            phase_list.extend(["p", "P"])
        if windows[0][0]["type"] == b"Swin":
            phase_list.extend(["s", "S"])
        if phase_list:
            min_dist = max(0, min(dist_all) - 10.0)
            max_dist = min(180, max(dist_all) + 10.0)
            dist_ttcurve = np.arange(min_dist, max_dist, 0.5)
            ttcurve = {}
            for phase_name in phase_list:
                ttcurve[phase_name] = []
            for dist in dist_ttcurve:
                arrivals = taup_model.get_travel_times(
                    source_depth_in_km=evdp,
                    distance_in_degree=dist,
                    phase_list=phase_list,
                )
                for arr in arrivals:
                    for phase_name in phase_list:
                        if arr.name == phase_name:
                            ttcurve[phase_name].append(
                                (arr.distance, arr.time, arr.ray_param)
                            )
            # sort (dist, ttime, rayp) points based on ray parameter
            for phase_name in phase_list:
                ttcurve[phase_name] = sorted(ttcurve[phase_name], key=lambda x: x[2])

        # ------ plot waveforms (one figure for each azimuthal bin of measurment windows)
        windows = sorted(
            windows, key=lambda x: (x[1]._v_attrs["azimuth"]) % 360
        )  # sort windows by azimuth
        nwin = len(windows)
        bin_idx0 = 0  # start index of windowes in the bin

        CRS = getattr(ccrs, config["plot"]["map_type"])
        map_params = config["plot"]["map_params"]
        projection = CRS(**map_params)
        if "map_extent" not in config["plot"]:
            min_lat = min(evla, min(stla_all)) - 2
            max_lat = max(evla, max(stla_all)) + 2
            min_lon = min(evlo, min(stlo_all)) - 2
            max_lon = max(evlo, max(stlo_all)) + 2
            map_extent = [min_lon, max_lon, min_lat, max_lat]
        else:
            map_extent = config["plot"]["map_extent"]

        min_weight = config["plot"]["min_weight"]
        min_SNR = min(config["misfit"]["window_weight"]["SNR"])

        with PdfPages(out_fig) as pdf:
            while bin_idx0 < nwin:
                winfo, g_sta = windows[bin_idx0]
                azmin = (g_sta._v_attrs["azimuth"]) % 360
                azmax = azmin + plot_azbin

                # get windows of one azimuthal bin
                windows_bin = []
                nwin_bin = 0
                while bin_idx0 < nwin:
                    winfo, g_sta = windows[bin_idx0]
                    bin_idx0 += 1 # index of next window
                    if not winfo["valid"]:
                        msg = f"invalid window: {winfo}, skip"
                        warnings.warn(msg)
                        continue
                    az = (g_sta._v_attrs["azimuth"]) % 360
                    dist_degree = g_sta._v_attrs["dist_degree"]
                    win_SNR = winfo["SNR"]
                    win_cc0 = winfo["cc0"]
                    win_ccmax = winfo["cc_max"]
                    # skip window not satisfying selection conditions
                    if plot_distlim.any():
                        if dist_degree < np.min(plot_distlim) or dist_degree > np.max(
                            plot_distlim
                        ):
                            continue
                    if plot_SNR and win_SNR < np.min(plot_SNR):
                        continue
                    if plot_CC0 and win_cc0 < np.min(plot_CC0):
                        continue
                    if plot_CCmax and win_ccmax < np.min(plot_CCmax):
                        continue
                    # put window into current bin or leave it for the next bin
                    if az <= azmax and (
                        not plot_max_ntrace_per_bin
                        or nwin_bin <= plot_max_ntrace_per_bin
                    ):
                        windows_bin.append((winfo, g_sta))
                        nwin_bin = nwin_bin + 1
                    else:
                        bin_idx0 -= 1
                        break

                # skip empty azimuthal_bin
                if nwin_bin == 0:
                    msg = f"No traces in the azimuthal bin [{bin_azmin}, {bin_azmax}], skip"
                    warnings.warn(msg)
                    continue

                bin_azmin = (windows_bin[0][1]._v_attrs["azimuth"]) % 360
                bin_azmax = (windows_bin[-1][1]._v_attrs["azimuth"]) % 360

                # ---- create figure
                print(
                    f"[INFO] plot azimuthal bin: {bin_azmin:05.1f} - {bin_azmax:05.1f}"
                )

                # fig = plt.figure(figsize=(8.5, 11)) # US letter
                fig = plt.figure(figsize=(8.27, 11.69))  # A4
                str_title = "{:s} ({:s} az:{:04.1f}~{:04.1f} dep:{:.1f})".format(
                    evnm, plot_window_id, bin_azmin, bin_azmax, evdp
                )
                if title_prefix:
                    str_title = "{:s} {:s}".format(title_prefix, str_title)
                fig.suptitle(str_title)
                # fig.text(
                #    0.5, 0.965, str_title, size="x-large", horizontalalignment="center"
                # )

                stla_bin = np.array(
                    [sta._v_attrs["latitude"] for _, sta in windows_bin]
                )
                stlo_bin = np.array(
                    [sta._v_attrs["longitude"] for _, sta in windows_bin]
                )
                weight_bin = np.array([win["weight"] for win, _ in windows_bin])
                ccdt_bin = np.array([win["cc_time_shift"] for win, _ in windows_bin])
                cc0_bin = np.array([win["cc0"] for win, _ in windows_bin])

                # event,station map with cc_dt
                ax_size = [0.35, 0.4]
                ax_origin = [0.02, 0.55]
                ax = fig.add_axes(ax_origin + ax_size, projection=projection)
                ax.set_extent(map_extent, crs=ccrs.PlateCarree())
                ax.coastlines(resolution="50m", color="black", linewidth=0.5)
                ax.add_feature(cfeature.BORDERS, linewidth=0.2, alpha=0.5)
                ax.gridlines(
                    draw_labels={"bottom": "x", "left": "y"},
                    formatter_kwargs={
                        "number_format": ".0f",
                        "degree_symbol": "",
                        "direction_label": False,
                    },
                    xlabel_style={"fontsize": 7},
                    ylabel_style={"fontsize": 7},
                    xpadding=1.5,
                    ypadding=1.5,
                    linewidth=0.2,
                )
                # ax.stock_img()
                ax.set_title(f"lat:{evla:.3f} lon:{evlo:.3f} dep:{evdp:.1f}")
                max_size = 6
                # sizes = max_size * weight_all**0.5
                # sizes[weight_all < min_weight] = max_size * min_weight**0.5
                im = ax.scatter(
                    stlo_all,
                    stla_all,
                    # s=sizes,
                    s=max_size,
                    # s=2,
                    c=ccdt_all,
                    marker="o",
                    cmap="seismic",
                    edgecolors="none",
                    # facecolors='none',
                    linewidth=0.1,
                    vmin=-10,
                    vmax=10,
                    # alpha=0.3,
                    transform=ccrs.Geodetic(),
                )
                # max_size = 6
                # sizes = max_size * weight_bin**0.5
                # sizes[weight_bin < min_weight] = max_size * min_weight**0.5
                # sizes = 3 * weight_bin
                # sizes[sizes < 1] = 1
                im = ax.scatter(
                    stlo_bin,
                    stla_bin,
                    # s=sizes,
                    s=max_size,
                    # s=5,
                    c=ccdt_bin,
                    marker="o",
                    cmap="seismic",
                    edgecolors="green",
                    linewidth=0.8,
                    vmin=-10,
                    vmax=10,
                    transform=ccrs.Geodetic(),
                )
                fig.colorbar(
                    im,
                    ax=ax,
                    location="bottom",
                    extend="both",
                    shrink=0.5,
                    fraction=0.1,
                    aspect=40,
                    label="cc_dt[sec]",
                )
                # plot focal mechanism
                sx, sy = projection.transform_point(evlo, evla, ccrs.Geodetic())
                # bb_width = 110000.0 * np.abs(max(stlo_all)-min(stlo_all)) * 0.1
                # bb_width = max(map_width, map_height) * 0.05
                x0, x1 = ax.get_xlim()
                bb_width = 0.05 * abs(x1 - x0)
                b = beach(
                    focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor="r"
                )
                ax.add_collection(b)

                # event,station map with cc0
                ax_size = [0.35, 0.4]
                ax_origin = [0.02, 0.05]
                ax = fig.add_axes(ax_origin + ax_size, projection=projection)
                ax.set_extent(map_extent, crs=ccrs.PlateCarree())
                ax.gridlines()
                ax.coastlines(linewidth=0.2)
                ax.stock_img()
                max_size = 6
                # sizes = max_size * weight_all**0.5
                # sizes[weight_all < min_weight] = max_size * min_weight**0.5
                im = ax.scatter(
                    stlo_all,
                    stla_all,
                    # s=sizes,
                    s=max_size,
                    # s=2,
                    c=cc0_all,
                    marker="o",
                    cmap="seismic",
                    edgecolors="none",
                    # facecolors='none',
                    linewidth=0.1,
                    vmin=0,
                    vmax=1,
                    # alpha=0.3,
                    transform=ccrs.Geodetic(),
                )
                # max_size = 6
                sizes = max_size * weight_bin**0.5
                sizes[weight_bin < min_weight] = max_size * min_weight**0.5
                im = ax.scatter(
                    stlo_bin,
                    stla_bin,
                    # s=sizes,
                    s=max_size,
                    # s=5,
                    c=cc0_bin,
                    marker="o",
                    cmap="seismic",
                    edgecolors="green",
                    linewidth=0.8,
                    vmin=0,
                    vmax=1,
                    transform=ccrs.Geodetic(),
                )
                fig.colorbar(
                    im,
                    ax=ax,
                    location="bottom",
                    extend="both",
                    shrink=0.5,
                    fraction=0.1,
                    aspect=40,
                    label="cc0",
                )
                # plot focal mechanism
                sx, sy = projection.transform_point(evlo, evla, ccrs.Geodetic())
                # bb_width = 110000.0 * np.abs(max(stlo_all)-min(stlo_all)) * 0.1
                # bb_width = max(map_width, map_height) * 0.05
                x0, x1 = ax.get_xlim()
                bb_width = 0.05 * abs(x1 - x0)
                b = beach(
                    focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor="r"
                )
                ax.add_collection(b)

                # create axis for seismograms
                ax_origin = [0.42, 0.06]
                ax_size = [0.4, 0.88]
                # ax_size = [0.3, 0.90]
                ax_1comp = fig.add_axes(ax_origin + ax_size)

                # -- xlim setting
                winb_bin = np.array(
                    [UTCDateTime(win["starttime"]) - evt0 for win, _ in windows_bin]
                )
                wine_bin = np.array(
                    [UTCDateTime(win["endtime"]) - evt0 for win, _ in windows_bin]
                )
                dist_bin = np.array(
                    [sta._v_attrs["dist_degree"] for _, sta in windows_bin]
                )
                # get time window relative to the regressed window central time
                plot_t0 = np.min(winb_bin - plot_rayp * dist_bin)
                plot_t1 = np.max(wine_bin - plot_rayp * dist_bin)
                plot_time = np.array([begin_time + plot_t0, end_time + plot_t1])

                # -- ylim setting
                y = dist_bin
                ny = len(y)
                plot_dy = 0.5 * (max(y) - min(y) + 1) / ny
                if plot_distlim.any():
                    plot_ymax = max(plot_distlim) + 2 * plot_dy
                    plot_ymin = min(plot_distlim) - 2 * plot_dy
                else:
                    plot_ymax = max(y) + 2 * plot_dy
                    plot_ymin = min(y) - 2 * plot_dy

                # -- plot traveltime curves
                for phase_name in phase_list:
                    # skip if no tt curves for this phase_names
                    if not ttcurve[phase_name]:
                        continue
                    # reduced time
                    phase_times = np.array(
                        [x[1] - plot_rayp * x[0] for x in ttcurve[phase_name]]
                    )
                    phase_distances = np.array([x[0] for x in ttcurve[phase_name]])
                    # skip if not in plot range
                    max_dist = np.max(phase_distances)
                    min_dist = np.min(phase_distances)
                    if max_dist < plot_ymin or min_dist > plot_ymax:
                        continue
                    ax_1comp.plot(phase_times, phase_distances, "b-", linewidth=0.1)
                    # ax_1comp.plot(phase_times, phase_distances, 'b.', markersize=0.5)
                    # label phase names
                    if max_dist < plot_ymax:
                        y_str = max_dist
                        x_str = max(phase_times[phase_distances == max_dist])
                    else:
                        y_str = plot_ymax
                        max_dist = max(phase_distances[phase_distances <= plot_ymax])
                        x_str = max(phase_times[phase_distances == max_dist])
                    ax_1comp.text(
                        x_str,
                        y_str,
                        phase_name,
                        verticalalignment="top",
                        horizontalalignment="center",
                        fontsize=11,
                        color="blue",
                    )

                # plot trace in each measurement window
                for win, g_sta in windows_bin:
                    attrs = g_sta._v_attrs
                    net = attrs["network"]
                    sta = attrs["station"]
                    stnm = f"{net}_{sta}"
                    first_arrtime = attrs["first_arrtime"]
                    dist_degree = attrs["dist_degree"]

                    print(f"plotting {stnm}")

                    obs = g_sta[obs_tag]
                    syn = g_sta[syn_tag]

                    sta_attrs = {
                        k: g_sta._v_attrs[k] for k in g_sta._v_attrs._f_list("user")
                    }
                    obs_attrs = {k: obs.attrs[k] for k in obs.attrs._f_list("user")}
                    syn_attrs = {k: syn.attrs[k] for k in syn.attrs._f_list("user")}
                    try:
                        # obs_ENZ, syn_ENZ, no_Hchan = self._extract_obs_syn_ENZ(
                        #     g_sta, obs_tag, syn_tag, event["tau"]
                        # )
                        obs_ENZ, syn_ENZ, no_Hchan = _extract_obs_syn_ENZ(
                            obs[:],
                            obs_attrs,
                            syn[:],
                            syn_attrs,
                            # sta_attrs, obs[:], obs_attrs, syn[:], syn_attrs
                        )
                        # debug
                        # np.savez(f'{net}_{sta}_old.npz', obs_ENZ=obs_ENZ, syn_ENZ=syn_ENZ, no_Hchan=no_Hchan)
                    except Exception as e:
                        msg = f"{stnm}: failed to get obs_ENZ,syn_ENZ ({e})"
                        warnings.warn(msg)
                        continue

                    obs_type = obs_attrs["type"]
                    syn_type = syn_attrs["type"]
                    if obs_type == "VEL" and syn_type == "DISP":
                        syn_to_vel = True
                    elif obs_type == "DISP" and syn_type == "DISP":
                        syn_to_vel = False
                    else:
                        msg = f"{stnm}: unknown obs/syn types ({obs_type}/{syn_type})"
                        raise ValueError(msg)
                    conv_stf = False
                    if syn_attrs["is_grn"]:
                        conv_stf = True
                        assert syn_attrs["origin_time"] == evt0

                    data_tb = obs.attrs["starttime"]
                    data_fs = obs.attrs["sampling_rate"]
                    data_dt = 1.0 / data_fs
                    data_nt = obs.attrs["npts"]
                    # data_te = data_tb + (data_nt - 1) * data_dt
                    # data_times = np.arange(data_nt, dtype=float) * data_dt
                    # noise_te = (first_arrtime - data_tb) - cfg_noise_before_first_arrival
                    # noise_idx1 = int(noise_te * data_fs)  # 0:idx1 as noise

                    # apply bandpass filter
                    butter_N = win["butter_N"]
                    butter_Wn = win["butter_Wn"]
                    # fft
                    npad = int(2 * data_fs / min(butter_Wn))
                    nfft = scipy.fft.next_fast_len(data_nt + npad)
                    freqs = np.fft.rfftfreq(nfft, d=data_dt)
                    # window filter
                    sos = scipy.signal.butter(
                        butter_N, butter_Wn, "bandpass", fs=data_fs, output="sos"
                    )
                    _, filter_h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
                    Fw = abs(filter_h)  # zero-phase filter response
                    # time derivative in frequency domain
                    Ft = 2j * np.pi * freqs

                    # obs = Fw * (Fd * d), obs_ENZ = Fd * d, d = disp. or vel.
                    obs_ENZ_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ, nfft), nfft)[
                        :, :data_nt
                    ]
                    # syn = Fw * (Fd * [Ft] * [S] * u), (syn_ENZ = Fd * u)
                    f_syn = np.fft.rfft(syn_ENZ, nfft)
                    if syn_to_vel:
                        f_syn *= Ft
                    if conv_stf:
                        f_src = stf_gauss_spectrum(freqs, evtau)
                        f_syn *= f_src
                    syn_ENZ_filt = np.fft.irfft(Fw * f_syn, nfft)[:, :data_nt]

                    # project to polarity defined by the window
                    cmpaz = win["cmpaz"]
                    cmpdip = win["cmpdip"]
                    sin_az = np.sin(np.deg2rad(cmpaz))
                    cos_az = np.cos(np.deg2rad(cmpaz))
                    sin_dip = np.sin(np.deg2rad(cmpdip))
                    cos_dip = np.cos(np.deg2rad(cmpdip))
                    cmp_vec = np.array(
                        [
                            cos_dip * sin_az,  # cos(E, comp)
                            cos_dip * cos_az,  # cos(N, comp)
                            -sin_dip,  # cos(Z, comp)
                        ]
                    )
                    obs_filt_proj = np.dot(cmp_vec, obs_ENZ_filt)
                    syn_filt_proj = np.dot(cmp_vec, syn_ENZ_filt)

                    # get plot time
                    reduced_time = dist_degree * plot_rayp
                    # time of first sample referred to centroid time
                    t0 = data_tb - evt0
                    # time of samples referred to centroid time
                    times = np.arange(data_nt) * data_dt + t0
                    # plot time window
                    plot_t0 = min(plot_time) + reduced_time
                    plot_t1 = max(plot_time) + reduced_time
                    plot_idx = (times > plot_t0) & (times < plot_t1)
                    # plot time (reduced time)
                    t_plot = times[plot_idx] - reduced_time

                    #  window begin/end
                    win_starttime = UTCDateTime(win["starttime"]) - evt0
                    win_endtime = UTCDateTime(win["endtime"]) - evt0
                    win_t0 = win_starttime - reduced_time
                    win_t1 = win_endtime - reduced_time
                    # win_idx = (times > win_starttime) & (times < win_endtime)
                    win_tshift = win["cc_time_shift"]

                    win_weight = win["weight"]
                    win_SNR = win["SNR"]

                    # plot seismograms
                    Amax_obs = win["obs_maxamp"]
                    Amax_syn = win["syn_maxamp"]

                    # # clip large amplitudes
                    # if plot_adj:
                    #     adj = station["adj"]
                    #     y = adj[plot_idx] / Amax_adj
                    #     idx = abs(y) > plot_clip + 1.0e-3
                    #     y[idx] = np.nan
                    #     ax_1comp.plot(
                    #         t_plot,
                    #         plot_flip * plot_dy * y + dist_degree,
                    #         "k-",
                    #         linewidth=0.5,
                    #     )

                    alpha = 1
                    # if win_weight < min_weight:
                    if win_SNR < min_SNR:
                        alpha = 0.1

                    y = obs_filt_proj[plot_idx] / Amax_obs
                    idx = abs(y) > plot_clip + 1.0e-3
                    y[idx] = np.nan
                    ax_1comp.plot(
                        t_plot,
                        plot_flip * plot_dy * y + dist_degree,
                        "k-",
                        linewidth=0.5,
                        alpha=alpha,
                    )

                    y = syn_filt_proj[plot_idx] / Amax_syn
                    idx = abs(y) > plot_clip + 1.0e-3
                    y[idx] = np.nan
                    ax_1comp.plot(
                        t_plot,
                        plot_flip * plot_dy * y + dist_degree,
                        "r-",
                        linewidth=0.5,
                        alpha=alpha,
                    )
                    if align_time:
                        ax_1comp.plot(
                            t_plot + win_tshift,
                            plot_flip * plot_dy * y + dist_degree,
                            "r--",
                            linewidth=0.5,
                            alpha=alpha,
                        )

                    # mark measure window range
                    ax_1comp.plot(win_t0, dist_degree, "b|", markersize=8, alpha=alpha)
                    ax_1comp.plot(win_t1, dist_degree, "b|", markersize=8, alpha=alpha)
                    ## annotate amplitude
                    #  ax.text(max(plot_time), dist_degree, '%.1e ' % (Amax_obs),
                    #      verticalalignment='bottom',
                    #      horizontalalignment='right',
                    #      fontsize=7, color='black')
                    #  ax.text(max(plot_time), dist_degree, '%.1e ' % (Amax_syn),
                    #      verticalalignment='top',
                    #      horizontalalignment='right',
                    #      fontsize=7, color='red')
                    ## annotate CC0
                    #  ax.text(max(plot_time), dist_degree, ' %.3f'%(window['cc']['CC0']),
                    #      verticalalignment='center', fontsize=7)
                    ## annotate window weight
                    # if i == 1:
                    #  ax.text(max(plot_time), dist_degree, ' %.1f' % (window['weight']),
                    #      verticalalignment='center', fontsize=7)
                    ##annotate station names
                    str_annot = "%s.%s(%.2f,%.1f,%.1f,%.1f)" % (
                        net,
                        sta,
                        win["cc0"],
                        win["cc_time_shift"],
                        win["SNR"],
                        win["weight"],
                    )
                    ax_1comp.text(
                        max(plot_time),
                        dist_degree,
                        str_annot,
                        verticalalignment="center",
                        fontsize=7,
                        alpha=alpha,
                    )
                    # ax_1comp.text(160, dist_degree, str_annot,
                    #    verticalalignment='center', fontsize=7)

                # -- set axes limits and lables, annotation
                ax_1comp.set_xlim(min(plot_time), max(plot_time))
                # ax_1comp.set_xlim(80,160)
                ax_1comp.set_ylim(plot_ymin, plot_ymax)
                ax_1comp.set_xlabel("t - {:.1f}*dist (s)".format(plot_rayp))
                ax_1comp.tick_params(axis="both", labelsize=10)
                # ylabel
                ax_1comp.set_ylabel("dist (deg)")
                # ax_1comp.invert_yaxis()

                # # -- save figures
                # if savefig:
                #     out_file = "%s/%s_az_%03d_%03d_%s.pdf" % (
                #         out_dir,
                #         evnm,  # event["id"],
                #         azmin,
                #         azmax,
                #         plot_window_id,
                #     )
                #     plt.savefig(out_file, format="pdf")
                # else:
                #     plt.show()
                pdf.savefig(fig)
                plt.close(fig)

        h5f.close()

    def output_adj(self, out_dir="adj"):

        h5f = pt.open_file(self.h5_path, "r")

        config = h5f.root._v_attrs["config"]
        adj_tag = config["adj"]["tag"]
        # band_code = config["adj"]["band_code"]
        suffix = config["adj"]["suffix"]

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])

        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        # TODO: parallelize by pool and imap_unordered
        for g_sta in g_wav:
            if adj_tag not in g_sta:
                continue
            adj_src = g_sta[adj_tag]

            tr = Trace()
            tr.stats.network = g_sta._v_attrs["network"]
            tr.stats.station = g_sta._v_attrs["station"]
            tr.stats.location = g_sta._v_attrs["location"]

            # time samples
            data_tb = adj_src.attrs["starttime"]
            data_fs = adj_src.attrs["sampling_rate"]
            data_dt = 1.0 / data_fs
            data_nt = adj_src.attrs["npts"]
            data_te = data_tb + (data_nt - 1) * data_dt

            solver_tb = adj_src.attrs["solver_tb"]
            solver_nt = adj_src.attrs["solver_nt"]
            solver_dt = adj_src.attrs["solver_dt"]
            solver_fs = 1.0 / solver_dt
            solver_te = solver_tb + (solver_nt - 1) * solver_dt

            npad_left = npad_right = 0
            if solver_tb < data_tb:
                npad_left = int((data_tb - solver_tb) * data_fs) + 1
            if solver_te > data_te:
                npad_right = int((solver_te - data_te) * data_fs) + 1

            solver_t0 = solver_tb - event_t0
            solver_times = (
                np.arange(solver_nt) * solver_dt + solver_t0
            )  # relative to origin time

            # write
            channels = adj_src.attrs["channels"]
            nchan = adj_src.shape[0]
            dat = np.zeros((solver_nt, 2), dtype=float)
            dat[:, 0] = solver_times
            for i in range(nchan):
                tr.data = np.zeros(data_nt + npad_left + npad_right, dtype=float)
                tr.data[npad_left : (npad_left + data_nt)] = adj_src[i, :]
                tr.stats.starttime = data_tb - npad_left * data_dt
                tr.stats.sampling_rate = data_fs
                cha = channels[i]["name"].decode()
                tr.stats.channel = cha

                out_file = os.path.join(out_dir, f"{tr.id}{suffix}")

                # tr.write(out_file + ".sac", format="sac")

                if solver_fs < data_fs:
                    msg = f"solver_fs ({solver_fs}) < data_fs ({data_fs})"
                    warnings.warn(msg)
                tr.interpolate(
                    starttime=solver_tb,
                    sampling_rate=solver_fs,
                    npts=solver_nt,
                    method="lanczos",
                    a=20,
                )

                # tr.write(out_file + ".1.sac", format="sac")

                # ascii format (needed by SEM)
                dat[:, 1] = tr.data
                np.savetxt(out_file, dat, fmt="%16.9e")

        h5f.close()

    def read_srcfrechet(
        self,
        srcfrechet_file,
    ):
        h5f = pt.open_file(self.h5_path, "r+")

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        # evnm = event['id'].decode()

        # read source gradient (src_frechet.000001)
        with open(srcfrechet_file, "r") as f:
            lines = [x for x in f.readlines() if not (x.startswith("#"))]

        lines = [x.split() for x in lines]

        # these outputs are not used
        # t0 = float(lines[0][0])
        # dchi_dt0 = float(lines[0][1])
        # tau = float(lines[1][0])
        # dchi_dtau = float(lines[1][1])

        xs = np.zeros(3)
        dchi_dxs = np.zeros(3)
        xs[0] = float(lines[2][0])
        dchi_dxs[0] = float(lines[2][1])
        xs[1] = float(lines[3][0])
        dchi_dxs[1] = float(lines[3][1])
        xs[2] = float(lines[4][0])
        dchi_dxs[2] = float(lines[4][1])

        mt = np.zeros((3, 3))
        dchi_dmt = np.zeros((3, 3))
        mt[0, 0] = float(lines[5][0])
        dchi_dmt[0, 0] = float(lines[5][1])
        mt[1, 1] = float(lines[6][0])
        dchi_dmt[1, 1] = float(lines[6][1])
        mt[2, 2] = float(lines[7][0])
        dchi_dmt[2, 2] = float(lines[7][1])

        mt[0, 1] = float(lines[8][0])
        dchi_dmt[0, 1] = float(lines[8][1])
        mt[1, 0] = mt[0, 1]
        dchi_dmt[1, 0] = dchi_dmt[0, 1]

        mt[0, 2] = float(lines[9][0])
        dchi_dmt[0, 2] = float(lines[9][1])
        mt[2, 0] = mt[0, 2]
        dchi_dmt[2, 0] = dchi_dmt[0, 2]

        mt[1, 2] = float(lines[10][0])
        dchi_dmt[1, 2] = float(lines[10][1])
        mt[2, 1] = mt[1, 2]
        dchi_dmt[2, 1] = dchi_dmt[1, 2]

        # correlation between mt and dchi_dmt
        # cc = np.sum(dchi_dmt * mt) / np.sum(mt**2) ** 0.5 / np.sum(dchi_dmt**2) ** 0.5
        # print("cc = ", cc)

        # check xs and mt are the same as in the event info
        assert np.allclose(xs, event["xs"])
        assert np.allclose(mt, event["mt"])

        event["grad_xs"] = dchi_dxs
        event["grad_mt"] = dchi_dmt
        tbl_src[0] = [event]
        tbl_src.flush()

        h5f.close()

    def make_perturbed_cmtsolution(
        self,
        out_dcmt_file="diff_CMTSOLUTION",
        out_cmtsolution_dxs="CMTSOLUTION_dxs",
        out_cmtsolution_dmt="CMTSOLUTION_dmt",
    ):
        """Make perturbed CMTSOLUTION based on source gradient"""
        h5f = pt.open_file(self.h5_path, "r+")

        config = h5f.root._v_attrs["config"]
        max_dxs_norm = config["source"]["max_dxs_norm"]
        max_dmt_ratio = config["source"]["max_dmt_ratio"]
        if max_dxs_norm <= 0.0:
            msg = f"max_dxs_norm ({max_dxs_norm}) must > 0"
            raise ValueError(msg)
        if max_dmt_ratio <= 0.0:
            msg = f"max_dmt_ratio ({max_dmt_ratio}) must > 0"
            raise ValueError(msg)

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        evhd = event["header"].decode()
        evnm = event["id"].decode()
        tau = event["tau"]
        xs = event["xs"]
        mt = event["mt"]
        grad_xs = event["grad_xs"]
        grad_xs /= np.max(np.abs(grad_xs))  # normalize in case of loss of precision
        grad_mt = event["grad_mt"]
        grad_mt /= np.max(
            np.abs(grad_mt)
        )  # normalize in case of loss of precision, e.g. for very small x, x*x would be 0

        # correlation between mt and dchi_dmt
        # cc = np.sum(dchi_dmt * mt) / np.sum(mt**2) ** 0.5 / np.sum(dchi_dmt**2) ** 0.5
        # print("cc = ", cc)

        # r0 = (np.sum(xs**2)) ** 0.5
        norm_grad_xs = (np.sum(grad_xs**2)) ** 0.5
        dxs_scaled = grad_xs
        if norm_grad_xs > 0:
            # dxs_scaled = max_dxs_ratio * r0 * grad_xs / norm_grad_xs
            dxs_scaled = max_dxs_norm * grad_xs / norm_grad_xs

        # make grad_mt orthogonal to mt, i.e. keep scalar moment unchanged
        grad_mt = grad_mt - mt * np.sum(grad_mt * mt) / np.sum(mt**2)
        m0 = (0.5 * np.sum(mt**2)) ** 0.5
        norm_grad_mt = (0.5 * np.sum(grad_mt**2)) ** 0.5
        dmt_scaled = grad_mt
        if norm_grad_mt > 0:
            dmt_scaled = max_dmt_ratio * m0 * grad_mt / norm_grad_mt
        # print(m0, norm_grad_mt, dmt_scaled)

        event["dxs"] = dxs_scaled
        event["dmt"] = dmt_scaled
        tbl_src[0] = [event]
        tbl_src.flush()

        h5f.close()

        # ====== write out differential CMTSOLUTION file
        with open(out_dcmt_file, "w") as fp:
            fp.write("%s\n" % evhd)
            fp.write("%-18s %s\n" % ("event name:", evnm))
            fp.write("%-18s %+15.8E\n" % ("dt0(s):", 0.0))
            fp.write("%-18s %+15.8E\n" % ("dtau(s):", 0.0))
            fp.write("%-18s %+.3f\n" % ("dx(m):", dxs_scaled[0]))
            fp.write("%-18s %+.3f\n" % ("dy(m):", dxs_scaled[1]))
            fp.write("%-18s %+.3f\n" % ("dz(m):", dxs_scaled[2]))
            fp.write("%-18s %+15.8E\n" % ("dMxx(N*m):", dmt_scaled[0, 0]))
            fp.write("%-18s %+15.8E\n" % ("dMyy(N*m):", dmt_scaled[1, 1]))
            fp.write("%-18s %+15.8E\n" % ("dMzz(N*m):", dmt_scaled[2, 2]))
            fp.write("%-18s %+15.8E\n" % ("dMxy(N*m):", dmt_scaled[0, 1]))
            fp.write("%-18s %+15.8E\n" % ("dMxz(N*m):", dmt_scaled[0, 2]))
            fp.write("%-18s %+15.8E\n" % ("dMyz(N*m):", dmt_scaled[1, 2]))

        # ====== write out perturbed CMTSOLUTION file
        xs_perturb = xs + dxs_scaled
        mt_perturb = mt + dmt_scaled
        # force mt_perturb to have the same scalar moment as mt
        mt_perturb = m0 * mt_perturb / (0.5 * np.sum(mt_perturb**2)) ** 0.5

        with open(out_cmtsolution_dxs, "w") as fp:
            fp.write("%s\n" % evhd)
            fp.write("%-18s %s_dxs\n" % ("event name:", evnm))
            fp.write("%-18s %+15.8E\n" % ("t0(s):", 0.0))
            fp.write("%-18s %+15.8E\n" % ("tau(s):", tau))
            fp.write("%-18s %+.3f\n" % ("x(m):", xs_perturb[0]))
            fp.write("%-18s %+.3f\n" % ("y(m):", xs_perturb[1]))
            fp.write("%-18s %+.3f\n" % ("z(m):", xs_perturb[2]))
            fp.write("%-18s %+15.8E\n" % ("Mxx(N*m):", mt[0, 0]))
            fp.write("%-18s %+15.8E\n" % ("Myy(N*m):", mt[1, 1]))
            fp.write("%-18s %+15.8E\n" % ("Mzz(N*m):", mt[2, 2]))
            fp.write("%-18s %+15.8E\n" % ("Mxy(N*m):", mt[0, 1]))
            fp.write("%-18s %+15.8E\n" % ("Mxz(N*m):", mt[0, 2]))
            fp.write("%-18s %+15.8E\n" % ("Myz(N*m):", mt[1, 2]))

        with open(out_cmtsolution_dmt, "w") as fp:
            fp.write("%s\n" % evhd)
            fp.write("%-18s %s_dmt\n" % ("event name:", evnm))
            fp.write("%-18s %+15.8E\n" % ("t0(s):", 0.0))
            fp.write("%-18s %+15.8E\n" % ("tau(s):", tau))
            fp.write("%-18s %+.3f\n" % ("x(m):", xs[0]))
            fp.write("%-18s %+.3f\n" % ("y(m):", xs[1]))
            fp.write("%-18s %+.3f\n" % ("z(m):", xs[2]))
            fp.write("%-18s %+15.8E\n" % ("Mxx(N*m):", mt_perturb[0, 0]))
            fp.write("%-18s %+15.8E\n" % ("Myy(N*m):", mt_perturb[1, 1]))
            fp.write("%-18s %+15.8E\n" % ("Mzz(N*m):", mt_perturb[2, 2]))
            fp.write("%-18s %+15.8E\n" % ("Mxy(N*m):", mt_perturb[0, 1]))
            fp.write("%-18s %+15.8E\n" % ("Mxz(N*m):", mt_perturb[0, 2]))
            fp.write("%-18s %+15.8E\n" % ("Myz(N*m):", mt_perturb[1, 2]))

    def update_source_dxs(self, dx, dy, dz):
        h5f = pt.open_file(self.h5_path, "r+")

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event["dxs"] = np.array([dx, dy, dz], dtype=float)
        tbl_src[0] = [event]
        tbl_src.flush()

        h5f.close()

    def linearized_search_cc(
        self,
        dm={"dt0": None, "dtau": None, "dxs": None, "dmt": None},
    ):
        """
        calculate normalized zero-lag cc for perturbed seismograms
        by linear combination of waveform derivatives of source parameters (t0,tau,xs,mt)

        dm={<model_name>:<model_step_sizes>, ...}

        return wcc_sum, weight_sum
        """
        # check model vectors in dm have the same length
        if (not dm) or len(dm) == 0:
            msg = "dm is empty!"
            raise ValueError(msg)
        # model_num = len(dm)

        # for model_name in dm:
        #     if model_name not in ["dt0", "dtau", "dxs", "dmt"]:
        #         msg = f"model_name[{model_name}] must be one of [t0, tau, xs, mt]"
        #         raise ValueError(msg)

        nstep = []
        for model_name in dm:
            if dm[model_name].ndim != 1:
                msg = f"dm[{model_name}] must be 1-D vector"
                raise ValueError(msg)
            nstep.append(np.size(dm[model_name]))

        if len(set(nstep)) != 1 or nstep[0] < 1:
            msg = "vectors in dm must have the same length"
            raise Exception(msg)
        nstep = nstep[0]

        h5f = pt.open_file(self.h5_path, "r")

        # check parameters
        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])

        if "dtau" in dm:
            tau = dm["dtau"] + event["tau"]
            if any(tau < 0):
                msg = f"too small dm['dtau'] ({dm['dtau']}) for event['tau']={event['tau']}"
                raise ValueError(msg)

        # syn_components = ["E", "N", "Z"]

        config = h5f.root._v_attrs["config"]
        obs_tag = config["data"]["tag"]
        syn_tag = config["syn"]["tag"]
        dsyn_grp = config["syn"]["dsyn_grp"]

        # sum of weighted normalized zero-lag cc at each model grid
        wcc_sum = np.zeros(nstep)
        # sum of all windows' weighting
        weight_sum = 0.0

        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        if "/window" not in h5f:
            msg = '"/window" not existing, run setup_window/measure_adj first!'
            raise KeyError(msg)
        tbl_win = h5f.get_node("/window")

        dsyn_tags = {}
        for model_name in dm:
            if model_name not in ["dt0", "dtau"]:
                dsyn_tags[model_name] = f"{dsyn_grp}/{model_name}"

        # print(dsyn_tags)

        for g_sta in g_wav:
            # print(f"cc_linearize for {g_sta._v_name}")

            attrs = g_sta._v_attrs
            net = attrs["network"]
            sta = attrs["station"]

            windows = tbl_win.read_where(
                f'(network == b"{net}") & (station == b"{sta}") & (valid == True) & (weight > 0)'
            )
            if len(windows) == 0:
                msg = f"No windows found for {g_sta._v_name}, SKIP"
                warnings.warn(msg)
                continue

            assert obs_tag in g_sta
            assert syn_tag in g_sta
            for _, tag in dsyn_tags.items():
                assert tag in g_sta

            obs = g_sta[obs_tag]
            syn = g_sta[syn_tag]

            obs_type = obs.attrs["type"]
            syn_type = syn.attrs["type"]
            if obs_type == "VEL" and syn_type == "DISP":
                syn_to_vel = True
            elif obs_type == "DISP" and syn_type == "DISP":
                syn_to_vel = False
            else:
                msg = (
                    f"unknown obs/syn types ({obs_type}/{syn_type}) in {g_sta._v_name}"
                )
                raise ValueError(msg)

            for _, dsyn_tag in dsyn_tags.items():
                assert g_sta[dsyn_tag].attrs["is_grn"] == g_sta[syn_tag].attrs["is_grn"]
                assert g_sta[dsyn_tag].attrs["type"] == g_sta[syn_tag].attrs["type"]

            try:
                dsyn_dm = {}
                obs_ENZ, no_Hchan = _get_obs_ENZ(g_sta, obs_tag)
                syn_ENZ = _get_syn_ENZ(g_sta, syn_tag)
                for m, tag in dsyn_tags.items():
                    dsyn_dm[m] = _get_syn_ENZ(g_sta, tag)
            except Exception as e:
                msg = f"failed to get {obs_tag},{syn_tag},{dsyn_tags} from {g_sta._v_name}, ({e})"
                warnings.warn(msg)
                continue

            if "dtau" in dm:
                assert g_sta[syn_tag].attrs["is_grn"] == True
                # for m, tag in dsyn_tags:
                #     assert g_sta[tag].attrs["is_grn"] == True

            conv_stf = False
            if g_sta[syn_tag].attrs["is_grn"]:
                conv_stf = True
                assert g_sta[syn_tag].attrs["origin_time"] == event_t0
                assert g_sta[dsyn_tag].attrs["origin_time"] == event_t0

            # dsyn_dm = {}
            # for model_name in dm:
            #     if model_name not in ["dt0", "dtau"]:
            #         dsyn_tag = f"{dsyn_grp}/{model_name}"
            #         try:
            #             dsyn_ENZ = self._get_syn_ENZ(g_sta, dsyn_tag)
            #         except Exception as e:
            #             msg = f"failed to get obs,syn_ENZ for {g_sta._v_name}, ({e})"
            #         dsyn_dm[model_name] = dsyn_ENZ

            # time samples
            data_tb = obs.attrs["starttime"]
            data_fs = obs.attrs["sampling_rate"]
            data_dt = 1.0 / data_fs
            data_nt = obs.attrs["npts"]
            data_times = np.arange(data_nt, dtype=float) * data_dt

            for win in windows:
                # print(win["id"])
                # skip bad windows
                if not win["valid"]:
                    msg = f"Window not measured for adj (WIN: {win}), SKIP"
                    warnings.warn(msg)
                    continue

                # window weight
                weight = win["weight"]

                # # skip window with zero weight
                # if np.isclose(weight, 0.0):
                #     msg = f"Window has near zero weight ({win['weight']}) (WIN: {win}), SKIP"
                #     warnings.warn(msg)
                #     continue

                weight_sum += weight

                # filter
                butter_N = win["butter_N"]
                butter_Wn = win["butter_Wn"]
                # fft
                npad = int(2 * data_fs / min(butter_Wn))
                nfft = scipy.fft.next_fast_len(data_nt + npad)
                freqs = np.fft.rfftfreq(nfft, d=data_dt)
                # window filter
                sos = scipy.signal.butter(
                    butter_N, butter_Wn, "bandpass", fs=data_fs, output="sos"
                )
                _, filter_h = scipy.signal.freqz_sos(sos, worN=freqs, fs=data_fs)
                Fw = abs(filter_h)  # zero-phase filter response
                # time derivative in frequency domain
                Ft = 2j * np.pi * freqs
                # time window
                win_tb = UTCDateTime(win["starttime"])
                win_te = UTCDateTime(win["endtime"])
                win_taper = win["taper"]
                # window function
                b = win_tb - data_tb
                e = win_te - data_tb
                win_len = e - b
                width = win_len * min(win_taper, 0.4)
                win_c = [b, b + width, e - width, e]
                win_func = cosine_sac_taper(data_times, win_c)
                # polarity projection
                proj_matrix = win["proj_matrix"]
                #
                obs_filt = np.fft.irfft(Fw * np.fft.rfft(obs_ENZ, nfft), nfft)[
                    :, :data_nt
                ]
                obs_filt_win = np.dot(proj_matrix, obs_filt) * win_func
                norm_obs = np.sqrt(np.sum(obs_filt_win**2))

                f_syn = np.fft.rfft(syn_ENZ, nfft)
                if syn_to_vel:
                    f_syn *= Ft
                f_dsyn_dm = {}
                for m, dsyn in dsyn_dm.items():
                    f_dsyn = np.fft.rfft(dsyn, nfft)
                    if syn_to_vel:
                        f_dsyn *= Ft
                    f_dsyn_dm[m] = f_dsyn

                # #DEBUG
                # f_syn0 = f_syn.copy()
                # if conv_stf:
                #     print(event["tau"])
                #     f_src = stf_gauss_spectrum(freqs, event["tau"])
                #     f_syn0 *= f_src
                # syn0 = np.fft.irfft(f_syn0, nfft)[:,:data_nt]
                # syn0_filt1 = np.fft.irfft(Fw * np.fft.rfft(syn0, nfft), nfft)[:,:data_nt]
                # syn0_filt = np.fft.irfft(Fw * f_syn0, nfft)[:,:data_nt]
                # syn0_filt_win = np.dot(proj_matrix, syn0_filt) * win_func
                # norm_syn0 = np.sqrt(np.sum(syn0_filt_win**2))
                # Nw = norm_obs * norm_syn0
                # cc0 = np.sum(obs_filt_win * syn0_filt_win) / Nw
                # np.savez('cc_linearized', obs_ENZ=obs_ENZ, syn_ENZ=syn_ENZ, syn0=syn0,
                #          syn0_filt=syn0_filt, syn0_filt1=syn0_filt1, win_func=win_func, obs_filt_win=obs_filt_win, syn_filt_win=syn0_filt_win)
                # print(cc0, win['cc0'])
                # assert np.isclose(cc0, win['cc0'])

                # search each model grid
                for istep in range(nstep):
                    f_syn1 = np.zeros_like(f_syn)
                    f_syn1 += f_syn
                    for m, f_dsyn in f_dsyn_dm.items():
                        f_syn1 += dm[m][istep] * f_dsyn

                    # perturbed source time function
                    dt0 = 0.0
                    if "dt0" in dm:
                        dt0 = dm["dt0"][istep]
                        f_tshift = np.exp(-2.0j * np.pi * freqs * dt0)
                        f_syn1 *= f_tshift

                    dtau = 0.0
                    if "dtau" in dm:
                        dtau = dm["dtau"][istep]
                    if conv_stf:
                        f_src = stf_gauss_spectrum(freqs, event["tau"] + dtau)
                        f_syn1 *= f_src

                    syn1_filt = np.fft.irfft(Fw * f_syn1, nfft)[:, :data_nt]
                    syn1_filt_win = np.dot(proj_matrix, syn1_filt) * win_func
                    norm_syn1 = np.sqrt(np.sum(syn1_filt_win**2))

                    # normalized cc between obs and perturbed syn
                    Nw = norm_obs * norm_syn1
                    cc = np.sum(obs_filt_win * syn1_filt_win) / Nw

                    # weighted cc
                    wcc_sum[istep] += weight * cc

        h5f.close()

        return wcc_sum, weight_sum

    def linearized_search_cc_multiprocess(
        self,
        dm={"dt0": None, "dtau": None, "dxs": None, "dmt": None},
        nproc=5,
    ):
        """
        calculate normalized zero-lag cc for perturbed seismograms
        by linear combination of waveform derivatives of source parameters (t0,tau,xs,mt)

        dm={<model_name>:<model_step_sizes>, ...}

        return wcc_sum, weight_sum
        """
        # check model vectors in dm have the same length
        if (not dm) or len(dm) == 0:
            msg = "dm is empty!"
            raise ValueError(msg)
        # model_num = len(dm)

        # for model_name in dm:
        #     if model_name not in ["dt0", "dtau", "dxs", "dmt"]:
        #         msg = f"model_name[{model_name}] must be one of [t0, tau, xs, mt]"
        #         raise ValueError(msg)

        nstep = []
        for model_name in dm:
            if dm[model_name].ndim != 1:
                msg = f"dm[{model_name}] must be 1-D vector"
                raise ValueError(msg)
            nstep.append(np.size(dm[model_name]))

        if len(set(nstep)) != 1 or nstep[0] < 1:
            msg = "vectors in dm must have the same length"
            raise Exception(msg)
        nstep = nstep[0]

        h5f = pt.open_file(self.h5_path, "r")

        # check parameters
        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        event_t0 = UTCDateTime(event["t0"])

        if "dtau" in dm:
            tau = dm["dtau"] + event["tau"]
            if any(tau < 0):
                msg = f"too small dm['dtau'] ({dm['dtau']}) for event['tau']={event['tau']}"
                raise ValueError(msg)

        if "/waveform" not in h5f:
            msg = '"/waveform" not existing, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        # get all station names
        stations = [g_sta._v_name for g_sta in g_wav]
        h5f.close()

        pool = multiprocessing.Pool(processes=nproc)
        inputs = [(self.h5_path, dm, stations[i::nproc]) for i in range(nproc)]
        results = pool.map(_linearized_search_cc_by_stations, inputs)

        wcc_sum = np.zeros(nstep)
        weight_sum = 0
        for wcc, weight in results:
            wcc_sum += wcc
            weight_sum += weight
        return wcc_sum, weight_sum

    def grid_search_source(self, out_txt, out_fig, max_niter=5, nproc=5):
        assert nproc >= 1

        range_shrink_ratio = 0.618

        dtau_range = 2
        dt0_range = 5
        dxs_range = 5
        dmt_range = 5

        dtau_opt = 0.0
        dt0_opt = 0.0
        dxs_opt = 0.0
        dmt_opt = 0.0

        # h5f = pt.open_file(self.h5_path, "r")
        with pt.open_file(self.h5_path, "r") as h5f:
            if "/source" not in h5f:
                msg = '"/source" not existing, run read_cmtsolution first!'
                raise KeyError(msg)
            tbl_src = h5f.get_node("/source")
            if tbl_src.nrows == 0:
                msg = "no source information"
                raise Exception(msg)
            event = tbl_src[0]
            tau0 = event["tau"]
        # h5f.close()

        f = open(out_txt, "w")

        niter = 0
        with PdfPages(out_fig) as pdf:
            while niter < max_niter:
                f.write("====== iter = {:02d}\n".format(niter))
                # print("DEBUG: iter = {:02d}\n".format(niter))

                # define search range
                dtau_min = dtau_opt - dtau_range
                if (tau0 + dtau_min) < 0.1:  # in case of negative tau
                    dtau_min = 0.1 - tau0
                dtau_max = dtau_opt + dtau_range
                dt0_min = dt0_opt - dt0_range
                dt0_max = dt0_opt + dt0_range
                dxs_min = dxs_opt - dxs_range
                dxs_max = dxs_opt + dxs_range
                dmt_min = dmt_opt - dmt_range
                dmt_max = dmt_opt + dmt_range

                f.write("search range for dtau: {:f} {:f}\n".format(dtau_min, dtau_max))
                f.write("search range for dt0:  {:f} {:f}\n".format(dt0_min, dt0_max))
                f.write("search range for dxs: {:f} {:f}\n".format(dxs_min, dxs_max))
                f.write("search range for dmt: {:f} {:f}\n".format(dmt_min, dmt_max))

                fig = plt.figure(figsize=(8.5, 11))
                fig.suptitle(f"iter#{niter:02d}")
                gs0 = fig.add_gridspec(2, 2)  # 2-D search between dt0 and dxs
                # gs1 = gs0[1,0].subgridspec(1, 2) # 1-D search of dtau and dmt

                # 2D grid search for dt0 and dxs
                dt0_1d = np.linspace(dt0_min, dt0_max, 5)
                dxs_1d = np.linspace(dxs_min, dxs_max, 5)
                dt0_2d, dxs_2d = np.meshgrid(dt0_1d, dxs_1d, indexing="ij")
                dm = {
                    "dt0": dt0_2d.reshape(dt0_2d.size),
                    "dxs": dxs_2d.reshape(dxs_2d.size),
                    "dmt": dmt_opt * np.ones(dt0_2d.size),
                    "dtau": dtau_opt * np.ones(dt0_2d.size),
                }
                # wcc_sum, weight_sum = self.linearized_search_cc(dm=dm)
                wcc_sum, weight_sum = self.linearized_search_cc_multiprocess(
                    dm=dm, nproc=nproc
                )
                # print(f'DEBUG: {wcc_sum=}, {weight_sum=}')

                # wcc = wcc_sum / weight_sum
                interp = scipy.interpolate.RectBivariateSpline(
                    dt0_1d, dxs_1d, wcc_sum.reshape(dt0_2d.shape)
                )
                dt0_1d_fine = np.linspace(np.min(dt0_1d), np.max(dt0_1d), 100)
                dxs_1d_fine = np.linspace(np.min(dxs_1d), np.max(dxs_1d), 100)
                wcc_sum_fine = interp(dt0_1d_fine, dxs_1d_fine)
                wcc_fine = wcc_sum_fine / weight_sum
                idx_max = np.argmax(wcc_sum_fine)
                ix, iy = np.unravel_index(idx_max, wcc_sum_fine.shape)
                wcc_sum_max = wcc_sum_fine[ix, iy]
                dt0_opt = dt0_1d_fine[ix]
                dxs_opt = dxs_1d_fine[iy]
                f.write("dt0_opt = {:f}\n".format(dt0_opt))
                f.write("dxs_opt = {:f}\n".format(dxs_opt))

                ax = fig.add_subplot(gs0[0, :])
                levels = np.linspace(np.min(wcc_fine), np.max(wcc_fine), 100)
                cs = ax.contour(dt0_1d_fine, dxs_1d_fine, wcc_fine.T, levels=levels)
                ax.clabel(cs, levels, inline=1, fmt="%.3f", fontsize=10)
                ax.plot(dt0_opt, dxs_opt, "ro", markersize=3, clip_on=False)
                ax.text(dt0_opt, dxs_opt, f"{dt0_opt:.1f},{dxs_opt:.1f}")
                ax.set_xlabel("dt0 (second)")
                ax.set_ylabel("dxs")
                ax.set_title(
                    f"max(wcc_sum)={wcc_sum_max:.2e}, weight_sum={weight_sum:.2e}"
                )

                # 1D grid search for dmt
                dmt_1d = np.linspace(dmt_min, dmt_max, 11)
                dm = {
                    "dmt": dmt_1d,
                    "dt0": dt0_opt * np.ones(dmt_1d.size),
                    "dxs": dxs_opt * np.ones(dmt_1d.size),
                    "dtau": dtau_opt * np.ones(dmt_1d.size),
                }
                # wcc_sum, weight_sum = self.linearized_search_cc(dm=dm)
                wcc_sum, weight_sum = self.linearized_search_cc_multiprocess(
                    dm=dm, nproc=nproc
                )
                wcc = wcc_sum / weight_sum
                interp = scipy.interpolate.interp1d(dmt_1d, wcc_sum, kind="cubic")
                dmt_1d_fine = np.linspace(dmt_min, dmt_max, 100)
                wcc_sum_fine = interp(dmt_1d_fine)
                wcc_fine = wcc_sum_fine / weight_sum
                idx_max = np.argmax(wcc_sum_fine)
                wcc_sum_max = wcc_sum_fine[idx_max]
                wcc_max = wcc_sum_max / weight_sum
                dmt_opt = dmt_1d_fine[idx_max]
                f.write("dmt_opt = {:f}\n".format(dmt_opt))

                ax1 = fig.add_subplot(gs0[1, 0])
                ax1.plot(dmt_1d, wcc, "ko")
                ax1.plot(dmt_1d_fine, wcc_fine, "k")
                ax1.plot(dmt_opt, wcc_max, "ro", markersize=5)
                ax1.text(dmt_opt, wcc_max, f"{dmt_opt:.1f}")
                ax1.set_xlabel("dmt")
                ax1.set_ylabel("wcc")
                ax1.set_title(f"max(wcc_sum)={wcc_sum_max:.2e}")

                # 1D grid search for dtau
                dtau_1d = np.linspace(dtau_min, dtau_max, 11)
                dm = {
                    "dtau": dtau_1d,
                    "dt0": dt0_opt * np.ones(dtau_1d.size),
                    "dxs": dxs_opt * np.ones(dtau_1d.size),
                    "dmt": dmt_opt * np.ones(dtau_1d.size),
                }
                # wcc_sum, weight_sum = self.linearized_search_cc(dm=dm)
                wcc_sum, weight_sum = self.linearized_search_cc_multiprocess(
                    dm=dm, nproc=nproc
                )
                wcc = wcc_sum / weight_sum
                interp = scipy.interpolate.interp1d(dtau_1d, wcc_sum, kind="cubic")
                dtau_1d_fine = np.linspace(dtau_min, dtau_max, 100)
                wcc_sum_fine = interp(dtau_1d_fine)
                wcc_fine = wcc_sum_fine / weight_sum
                idx_max = np.argmax(wcc_sum_fine)
                wcc_sum_max = wcc_sum_fine[idx_max]
                wcc_max = wcc_sum_max / weight_sum
                dtau_opt = dtau_1d_fine[idx_max]
                f.write("dtau_opt = {:f}\n".format(dtau_opt))

                ax2 = fig.add_subplot(gs0[1, 1], sharey=ax1)
                ax2.plot(dtau_1d, wcc, "ko")
                ax2.plot(dtau_1d_fine, wcc_fine, "k")
                ax2.plot(dtau_opt, wcc_max, "ro", markersize=5)
                ax2.text(dtau_opt, wcc_max, f"{dtau_opt:.1f}")
                ax2.set_xlabel("dtau")
                ax2.tick_params(labelleft=False)
                ax2.set_title(f"max(wcc_sum)={wcc_sum_max:.2e}")
                # ax2.set_ylabel("wcc(interp)")

                f.write(f"weight_sum = {weight_sum:.5e}\n")
                f.write(f"max(wcc_sum) = {wcc_sum_max:.5e}\n")
                f.write(f"max(wcc_sum)/weight_sum = {wcc_max:.5e}\n")

                pdf.savefig()
                plt.close()

                # shrink search range
                dtau_range = range_shrink_ratio * dtau_range
                dt0_range = range_shrink_ratio * dt0_range
                dxs_range = range_shrink_ratio * dxs_range
                dmt_range = range_shrink_ratio * dmt_range

                niter += 1

        # ====== close file
        f.close()

        # save final step lengths
        with pt.open_file(self.h5_path, "r+") as h5f:
            if "/source" not in h5f:
                msg = '"/source" not existing, run read_cmtsolution first!'
                raise KeyError(msg)
            tbl_src = h5f.get_node("/source")
            if tbl_src.nrows == 0:
                msg = "no source information"
                raise Exception(msg)
            event = tbl_src[0]
            event["step_t0"] = dt0_opt
            event["step_tau"] = dtau_opt
            event["step_xs"] = dxs_opt
            event["step_mt"] = dmt_opt
            tbl_src[0] = [event]
            tbl_src.flush()

    def make_updated_cmtsolution(
        self, cmtfile="CMTSOLUTION.updated", dt0=0.0, dtau=0.0, dxs=0.0, dmt=0.0
    ):
        h5f = pt.open_file(self.h5_path, "r")

        if "/source" not in h5f:
            msg = '"/source" not existing, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]
        evhd = event["header"].decode()
        evnm = event["id"].decode()
        t0 = UTCDateTime(event["t0"])
        tau = event["tau"]
        xs = event["xs"]
        mt = event["mt"]
        dxs_vec = event["dxs"]
        dmt_vec = event["dmt"]

        # update origin time in header line
        t1 = t0 + dt0
        header = evhd.split()
        header[1] = "{:04d}".format(t1.year)
        header[2] = "{:02d}".format(t1.month)
        header[3] = "{:02d}".format(t1.day)
        header[4] = "{:02d}".format(t1.hour)
        header[5] = "{:02d}".format(t1.minute)
        minisecond = round(t1.microsecond * 1e-3)  # precise to minisecond
        header[6] = "{:02d}.{:03d}".format(t1.second, minisecond)
        header_line = " ".join(header)
        # update source duration tau
        tau1 = tau + dtau
        # update source location
        xs1 = xs + dxs * dxs_vec
        # update moment tensor
        mt1 = mt + dmt * dmt_vec

        # force updated mt to have the same scalar moment as mt since absolute amplitude is ignored
        m0 = (0.5 * np.sum(mt**2)) ** 0.5
        mt1 = m0 * mt1 / (0.5 * np.sum(mt1**2)) ** 0.5

        # write updated cmt file
        with open(cmtfile, "w") as fp:
            fp.write("%s\n" % (header_line))
            fp.write("%-18s %s\n" % ("event name:", evnm))
            fp.write("%-18s %+15.8E\n" % ("t0(s):", 0.0))
            fp.write("%-18s %+15.8E\n" % ("tau(s):", tau1))
            fp.write("%-18s %+.3f\n" % ("x(m):", xs1[0]))
            fp.write("%-18s %+.3f\n" % ("y(m):", xs1[1]))
            fp.write("%-18s %+.3f\n" % ("z(m):", xs1[2]))
            fp.write("%-18s %+15.8E\n" % ("Mxx(N*m):", mt1[0, 0]))
            fp.write("%-18s %+15.8E\n" % ("Myy(N*m):", mt1[1, 1]))
            fp.write("%-18s %+15.8E\n" % ("Mzz(N*m):", mt1[2, 2]))
            fp.write("%-18s %+15.8E\n" % ("Mxy(N*m):", mt1[0, 1]))
            fp.write("%-18s %+15.8E\n" % ("Mxz(N*m):", mt1[0, 2]))
            fp.write("%-18s %+15.8E\n" % ("Myz(N*m):", mt1[1, 2]))

        h5f.close()

    def get_window_ids(self):
        with pt.open_file(self.h5_path, "r") as h5f:
            win_tbl = h5f.root["window"]
            win_ids = sorted([w.decode() for w in set(win_tbl.read(field="id"))])
        return win_ids

    def plot_misfit(self, out_fig, title_prefix=None):
        """Plot misfit for a certain event and window_id"""

        h5f = pt.open_file(self.h5_path, "r")

        # config
        config = h5f.root._v_attrs["config"]

        # event info
        if "/source" not in h5f:
            msg = '"/source" does not exist, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]

        if "/waveform" not in h5f:
            msg = '"/waveform" does not exist, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        if "/window" not in h5f:
            msg = '"/waveform" does not exist, run setup_windows first!'
            raise KeyError(msg)
        tbl_win = h5f.get_node("/window")

        # event info
        evt0 = UTCDateTime(event["t0"])
        evtau = event["tau"]
        evla = event["latitude"]
        evlo = event["longitude"]
        evdp = event["depth"]
        evnm = event["id"].decode()
        # evdp has to be >=0 otherwise taup would crash
        if evdp < 0.0:
            evdp = 0.0
        mt = event["mt_rtp"]
        Mrr = mt[0][0]
        Mtt = mt[1][1]
        Mpp = mt[2][2]
        Mrt = mt[0][1]
        Mrp = mt[0][2]
        Mtp = mt[1][2]
        focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]

        # get all windows
        data = [
            (win, g_wav[stnm])
            for win in tbl_win.read()
            if (stnm := f'{win["network"].decode()}_{win["station"].decode()}') in g_wav
        ]
        if not data:
            warnings.warn("No data to plot!")
            return

        # map configuration
        # parallels = np.arange(-90.0, 90, 10.0)
        # meridians = np.arange(0.0, 360, 10.0)

        stla_all = np.array([sta._v_attrs["latitude"] for _, sta in data])
        stlo_all = np.array([sta._v_attrs["longitude"] for _, sta in data])

        # map_latc = kwargs.get("center_lat", np.median(stla_all))
        # map_lonc = kwargs.get("center_lon", np.median(stlo_all))
        # min_lon = kwargs.get("min_lon", min(min(stlo_all), evlo))
        # max_lon = kwargs.get("max_lon", max(max(stlo_all), evlo))
        # min_lat = kwargs.get("min_lat", min(min(stla_all), evla))
        # max_lat = kwargs.get("max_lat", max(max(stla_all), evla))

        # map_lat0, map_lon0, map_width, map_height = centerMap(
        #     [*stla_all, evla], [*stlo_all, evlo], 1.1
        # )

        window_ids = list(set([win["id"].decode() for win, _ in data]))

        # projection = ccrs.TransverseMercator(central_latitude=map_latc, central_longitude=map_lonc)
        # projection = ccrs.Mercator(central_longitude=map_lonc, min_latitude=min_lat, max_latitude=max_lat,
        #                            latitude_true_scale=map_latc)
        # map_config['type']
        CRS = getattr(ccrs, config["plot"]["map_type"])
        map_params = config["plot"]["map_params"]
        projection = CRS(**map_params)

        with PdfPages(out_fig) as pdf:
            for window_id in window_ids:
                data_sel = [
                    (win, sta) for win, sta in data if win["id"].decode() == window_id
                ]
                stlats = np.array([sta._v_attrs["latitude"] for _, sta in data_sel])
                stlons = np.array([sta._v_attrs["longitude"] for _, sta in data_sel])
                ccdt_list = [win["cc_time_shift"] for win, _ in data_sel]
                cc0_list = [win["cc0"] for win, _ in data_sel]
                weight_list = [win["weight"] for win, _ in data_sel]

                # fig = plt.figure(figsize=(8.5, 11))
                fig = plt.figure(figsize=(8.27, 11.69))  # A4
                str_title = "{:s} ({:s} dep:{:.1f})".format(evnm, window_id, evdp)
                if title_prefix:
                    str_title = "{:s} {:s}".format(title_prefix, str_title)
                fig.suptitle(str_title)
                # fig.text(
                #     0.5, 0.965, str_title, size="x-large", horizontalalignment="center"
                # )

                # plots of cc_dt
                ax = fig.add_subplot(211, projection=projection)
                # ax.set_extent((min_lon, max_lon, min_lat, max_lat), crs=ccrs.Geodetic())

                # Draw standard features
                ax.gridlines()
                ax.coastlines()
                ax.stock_img()

                sizes = 5 * np.array(weight_list)
                sizes[sizes < 1] = 1
                im = ax.scatter(
                    stlons,
                    stlats,
                    s=sizes,
                    c=ccdt_list,
                    marker="o",
                    cmap="seismic",
                    edgecolor="gray",
                    linewidth=0.1,
                    vmin=-10,
                    vmax=10,
                    transform=ccrs.Geodetic(),
                )
                fig.colorbar(
                    im,
                    ax=ax,
                    location="right",
                    extend="both",
                    shrink=0.5,
                    fraction=0.1,
                    aspect=40,
                    label="cc_dt",
                )

                # plot focal mechanism
                sx, sy = projection.transform_point(evlo, evla, ccrs.Geodetic())
                # bb_width = 110000.0 * np.abs(max(stlo_all)-min(stlo_all)) * 0.1
                # bb_width = max(map_width, map_height) * 0.05
                x0, x1 = ax.get_xlim()
                bb_width = 0.05 * abs(x1 - x0)
                b = beach(
                    focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor="r"
                )
                ax.add_collection(b)

                # plots of cc0
                ax = fig.add_subplot(212, projection=projection)
                # ax.set_extent((min_lon, max_lon, min_lat, max_lat), crs=ccrs.Geodetic())

                # Draw standard features
                ax.gridlines()
                ax.coastlines()
                ax.stock_img()

                sizes = 5 * np.array(weight_list)
                sizes[sizes < 1] = 1
                im = ax.scatter(
                    stlons,
                    stlats,
                    s=sizes,
                    c=cc0_list,
                    marker="o",
                    cmap="seismic",
                    edgecolor="gray",
                    linewidth=0.1,
                    vmin=0,
                    vmax=1,
                    transform=ccrs.Geodetic(),
                )
                fig.colorbar(
                    im,
                    ax=ax,
                    location="right",
                    extend="both",
                    shrink=0.5,
                    fraction=0.1,
                    aspect=40,
                    label="cc0",
                )

                # plot focal mechanism
                sx, sy = projection.transform_point(evlo, evla, ccrs.Geodetic())
                x0, x1 = ax.get_xlim()
                bb_width = 0.05 * abs(x1 - x0)
                b = beach(
                    focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor="r"
                )
                ax.add_collection(b)

                pdf.savefig(fig)
                plt.close(fig)

        h5f.close()

    def plot_histogram(self, out_fig, title_prefix=None, min_SNR=5.0, max_dt=20.0):
        h5f = pt.open_file(self.h5_path, "r")

        # config
        config = h5f.root._v_attrs["config"]
        # min_weight = config['plot']['min_weight']
        nbins = config["plot"]["nbins"]
        if "cc_tshift" in config["misfit"]["window_weight"]:
            max_dt = max(config["misfit"]["window_weight"]["cc_tshift"])
        if "SNR" in config["misfit"]["window_weight"]:
            min_SNR = min(config["misfit"]["window_weight"]["SNR"])

        # event info
        if "/source" not in h5f:
            msg = '"/source" does not exist, run read_cmtsolution first!'
            raise KeyError(msg)
        tbl_src = h5f.get_node("/source")
        if tbl_src.nrows == 0:
            msg = "no source information"
            raise Exception(msg)
        event = tbl_src[0]

        if "/waveform" not in h5f:
            msg = '"/waveform" does not exist, run read_data_h5 first!'
            raise KeyError(msg)
        g_wav = h5f.get_node("/waveform")

        if "/window" not in h5f:
            msg = '"/waveform" does not exist, run setup_windows first!'
            raise KeyError(msg)
        tbl_win = h5f.get_node("/window")

        # event info
        evnm = event["id"].decode()
        evdp = event["depth"]

        # get all windows
        data = [
            (win, g_wav[stnm])
            for win in tbl_win.read()
            if (stnm := f'{win["network"].decode()}_{win["station"].decode()}') in g_wav
        ]
        if not data:
            warnings.warn("No data to plot!")
            return

        window_ids = sorted(set([win["id"].decode() for win, _ in data]))

        with PdfPages(out_fig) as pdf:
            snr_list = np.array([win["SNR"] for win, _ in data])
            ccdt_list = np.array([win["cc_time_shift"] for win, _ in data])
            cc0_list = np.array([win["cc0"] for win, _ in data])
            weight_list = np.array([win["weight"] for win, _ in data])
            wcc_sum = sum(cc0_list * weight_list)

            fig = plt.figure(figsize=(10, 5))
            str_title = "{:s} (dep:{:.1f}): all (sum(weight*cc0)={:.1f})".format(
                evnm, evdp, wcc_sum
            )
            if title_prefix:
                str_title = "{:s} {:s}".format(title_prefix, str_title)
            fig.suptitle(str_title)
            ax = fig.add_subplot(121)
            dt_used = ccdt_list[abs(ccdt_list) <= max_dt]
            mean_dt = np.mean(dt_used)
            std_dt = np.std(dt_used)
            label = f"mean(std) = {mean_dt:.3f}({std_dt:.1f})"
            ax.hist(dt_used, nbins, label=label)
            ax.set_xlabel("dt_cc [obs-syn] (second)")
            ax.set_xlim([-max_dt, max_dt])
            ax.legend()
            ax = fig.add_subplot(122)
            ax.hist(cc0_list[snr_list >= min_SNR], nbins)
            ax.set_xlabel("cc0")
            ax.set_xlim([0, 1])
            pdf.savefig(fig)
            plt.close(fig)

            for window_id in window_ids:
                data_sel = [
                    (win, sta) for win, sta in data if win["id"].decode() == window_id
                ]
                snr_list = np.array([win["SNR"] for win, _ in data_sel])
                ccdt_list = np.array([win["cc_time_shift"] for win, _ in data_sel])
                cc0_list = np.array([win["cc0"] for win, _ in data_sel])
                weight_list = np.array([win["weight"] for win, _ in data_sel])
                wcc_sum = sum(cc0_list * weight_list)

                fig = plt.figure(figsize=(10, 5))
                str_title = "{:s} (dep:{:.1f}): {:s} (sum(weight*cc0)={:.1f})".format(
                    evnm, evdp, window_id, wcc_sum
                )
                if title_prefix:
                    str_title = "{:s} {:s}".format(title_prefix, str_title)
                fig.suptitle(str_title)
                ax = fig.add_subplot(121)
                dt_used = ccdt_list[abs(ccdt_list) <= max_dt]
                mean_dt = np.mean(dt_used)
                std_dt = np.std(dt_used)
                label = f"mean(std) = {mean_dt:.3f}({std_dt:.1f})"
                ax.hist(dt_used, nbins, label=label)
                ax.set_xlabel("dt_cc [obs-syn] (second)")
                ax.set_xlim([-max_dt, max_dt])
                ax.legend()
                # ax.set_title(f"mean(std) = {mean_dt}({std_dt})")
                ax = fig.add_subplot(122)
                ax.hist(cc0_list[snr_list >= min_SNR], nbins)
                ax.set_xlabel("cc0")
                ax.set_xlim([0, 1])
                pdf.savefig(fig)
                plt.close(fig)

        h5f.close()

    def grid_search_structure(self, dm={"dm": np.linspace(0, 2, 20)}, nproc=5):
        """
        dm = {'tag1': np.linspace(b1,e1,n1),
              'tag2': np.linspace(b2,e2,n2),
              ...
              }
        """
        assert nproc >= 1

        # for storing data, /waveforms/NET_STA/DATA_DISP[nchan,nt]
        h5_atom = pt.Atom.from_dtype(np.dtype(np.float64))
        h5_filters = pt.Filters(complevel=3, complib="zlib")

        # with pt.open_file(self.h5_path, "r") as h5f:
        #     config = h5f.root._v_attrs["config"]
        #     grid_search = config['structure']['grid_search']

        # sampling points along each search direction
        # dm = {tag: np.linspace(*vals) for tag, vals in grid_search.items()}

        # create grid of search points expanded by all dm's
        grids = np.meshgrid(*dm.values(), indexing="ij")
        grid_shape = grids[0].shape
        dm_grid = {tag: grids[i].flatten() for i, tag in enumerate(dm)}

        # wcc_sum, weight_sum = self.linearized_search_cc(dm=dm)
        wcc_sum, weight_sum = self.linearized_search_cc_multiprocess(
            dm=dm_grid, nproc=nproc
        )
        wcc_sum = wcc_sum.reshape(grid_shape)

        # save dm, wcc_sum, weight_sum
        with pt.open_file(self.h5_path, "r+") as h5f:
            if "/grid_search/structure" not in h5f:
                g_search = h5f.create_group(
                    "/grid_search", "structure", createparents=True
                )
            else:
                g_search = h5f.get_node("/grid_search/structure")

            for tag in dm:
                if tag in g_search:
                    h5f.remove_node(g_search, tag)
                shape = (len(dm[tag]),)
                ca = h5f.create_carray(
                    g_search, tag, h5_atom, shape, filters=h5_filters
                )
                ca[:] = np.array(dm[tag], dtype=np.float64)

            if "wcc_sum" in g_search:
                h5f.remove_node(g_search, "wcc_sum")
            ca = h5f.create_carray(
                g_search, "wcc_sum", h5_atom, grid_shape, filters=h5_filters
            )
            ca[:] = np.array(wcc_sum, dtype=np.float64)
            ca.attrs["weight_sum"] = weight_sum
            ca.attrs["dims"] = list(
                dm.keys()
            )  # corresponding dm name of each dimension of wcc_sum


#
#     def output_misfit(self, out_file):
#         """
#         Output misfit measurements
#
#         """
#         event = self.data["event"]
#         station_dict = self.data["station"]
#
#         fp = open(out_file, "w")
#         fp.write("#station window weight CC0 CCmax dt_cc SNR AR0 ARmax\n")
#
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected stations
#             if station["stat"]["code"] < 1:
#                 continue
#
#             window_dict = station["window"]
#             for window_id in window_dict:
#                 window = window_dict[window_id]
#                 # skip rejected/un-measured windows
#                 if window["stat"]["code"] < 1:
#                     continue
#
#                 # skip bad windows
#                 # if window['weight'] < 1.0e-3:
#                 #  continue
#
#                 cc = window["cc"]
#                 quality = window["quality"]
#
#                 fp.write(
#                     "{:15s} {:15s} {:7.5f}  {:7.5f}  {:7.5f}  {:+7.3f}"
#                     "  {:7.3f}  {:12.5e}  {:12.5e}\n".format(
#                         station_id,
#                         window_id,
#                         window["weight"],
#                         cc["CC0"],
#                         cc["CCmax"],
#                         cc["cc_tshift"],
#                         quality["SNR"],
#                         cc["AR0"],
#                         cc["ARmax"],
#                     )
#                 )
#
#         fp.close()
#
#     #
#     # ======================================================
#     #
#
#     def plot_histogram(
#         self,
#         max_dt=10,
#         nbins=100,
#         outfig="cc_tshift_histogram.pdf",
#     ):
#         """
#         Plot histograms of measured cc time shifts for each windows
#
#         Parameters
#         ---------
#         (string) window_type: type of windows to plot
#             only surface or body are supported
#
#         (float) max_dt: maximum abs(dt) to use
#         """
#         event = self.data["event"]
#         station_dict = self.data["station"]
#
#         # get measurement results
#         data = {}
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected stations
#             if station["stat"]["code"] < 1:
#                 continue
#             window_dict = station["window"]
#             # loop all windows
#             for window_id in window_dict:
#                 window = window_dict[window_id]
#                 # skip rejected/un-measured windows
#                 if window["stat"]["code"] < 1:
#                     continue
#                 # skip bad windows
#                 if window["weight"] < 1.0e-3:
#                     continue
#                 # append info to list
#                 if window_id not in data:
#                     data[window_id] = []
#                 cc_dict = window["cc"]
#                 data[window_id].append(cc_dict["cc_tshift"])
#
#         # ------ plot
#         with PdfPages(outfig) as pdf:
#             # individual windows
#             for window_id in data:
#                 dt_cc = np.array(data[window_id])
#                 idx = np.abs(dt_cc) <= max_dt
#                 plt.hist(dt_cc[idx], nbins, histtype="step")
#                 title_str = "%s %s %f$\pm$%f" % (
#                     event["id"],
#                     window_id,
#                     np.mean(dt_cc[idx]),
#                     np.std(dt_cc[idx]),
#                 )
#                 plt.title(title_str)
#                 plt.xlabel("dt_cc [obs-syn] (second)")
#                 plt.ylabel("Window number")
#                 pdf.savefig()
#                 plt.close()
#             # body-wave windows
#             dt_cc = []
#             for window_id in data:
#                 if "surface_" not in window_id:
#                     dt_cc += data[window_id]
#             dt_cc = np.array(dt_cc)
#             idx = np.abs(dt_cc) <= max_dt
#             plt.hist(dt_cc[idx], nbins, histtype="step")
#             title_str = "%s body-wave %f$\pm$%f" % (
#                 event["id"],
#                 np.mean(dt_cc[idx]),
#                 np.std(dt_cc[idx]),
#             )
#             plt.title(title_str)
#             plt.xlabel("dt_cc [obs-syn] (second)")
#             plt.ylabel("Window number")
#             pdf.savefig()
#             plt.close()
#             # surface_RZ
#             dt_cc = []
#             for window_id in data:
#                 if ("surface_R" in window_id) or ("surface_Z" in window_id):
#                     dt_cc += data[window_id]
#             dt_cc = np.array(dt_cc)
#             idx = np.abs(dt_cc) <= max_dt
#             plt.hist(dt_cc[idx], nbins, histtype="step")
#             title_str = "%s surface_RZ %f$\pm$%f" % (
#                 event["id"],
#                 np.mean(dt_cc[idx]),
#                 np.std(dt_cc[idx]),
#             )
#             plt.title(title_str)
#             plt.xlabel("dt_cc [obs-syn] (second)")
#             plt.ylabel("Window number")
#             pdf.savefig()
#             plt.close()
#
#     #
#     # ======================================================
#     #
#
#     def plot_histogram_combine(
#         self,
#         window_type="surface",
#         max_dt=10,
#         nbins=100,
#         out_file="./histogram_dt_cc.pdf",
#     ):
#         """
#         Plot histograms of measured cc time shifts for combined window type
#
#         Parameters
#         ---------
#         (string) window_type: type of windows to plot
#             only surface or body are supported
#
#         (float) max_dt: maximum abs(dt) to use
#
#         """
#         event = self.data["event"]
#         station_dict = self.data["station"]
#
#         # get measurement results
#         weight = []
#         cc0 = []
#         dt_cc = []
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected stations
#             if station["stat"]["code"] < 1:
#                 continue
#             window_dict = station["window"]
#             # loop all windows
#             for window_id in window_dict:
#                 window = window_dict[window_id]
#                 # skip rejected/un-measured windows
#                 if window["stat"]["code"] < 1:
#                     continue
#                 # check window type
#                 if window_type == "surface":
#                     if "surface" not in window_id:
#                         continue
#                 elif window_type == "body":
#                     if "surface" in window_id:
#                         continue
#                 else:
#                     raise Exception("window_type should be either surface or body !")
#                 # skip bad windows
#                 if window["weight"] < 1.0e-3:
#                     continue
#                 # append info to list
#                 weight = window["weight"]
#                 cc_dict = window["cc"]
#                 cc0.append(cc_dict["CC0"])
#                 dt_cc.append(cc_dict["cc_tshift"])
#
#         # plot
#         weight = np.array(weight)
#         cc0 = np.array(cc0)
#         dt_cc = np.array(dt_cc)
#
#         idx = np.abs(dt_cc) <= max_dt
#
#         dt_mean = np.mean(dt_cc)
#         dt_std = np.std(dt_cc)
#
#         plt.hist(dt_cc[idx], nbins=nbins)
#         title_str = "%s %s_wave mean(std)=%f(%f)" % (
#             event["id"],
#             window_type,
#             dt_mean,
#             dt_std,
#         )
#         plt.title(title_str)
#         plt.xlabel("dt_cc [obs-syn] (second)")
#         plt.ylabel("Event number")
#         plt.savefig(out_file, format="pdf")
#
#     def waveform_der_dmodel(
#         self,
#         syn_dir="output_perturb/sac",
#         syn_band_code="MX",
#         syn_suffix=".sem.sac",
#         model_name="vsv",
#         sac_dir=None,
#     ):
#         """
#         Get partial derivative seismograms from synthetics of perturbed model
#         along the model update direction
#
#         Notes
#         -----
#         use finite difference to get waveform derivative
#
#         must keep all source parameter the same as used in the forward simulation.
#
#         """
#         syn_orientation_codes = ["E", "N", "Z"]
#
#         event = self.data["event"]
#         t0 = event["t0"]
#
#         station_dict = self.data["station"]
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected statations
#             if station["stat"]["code"] < 1:
#                 continue
#
#             # ------ time samples
#             waveform = station["waveform"]
#             time_sample = waveform["time_sample"]
#             starttime = time_sample["starttime"]
#             dt = time_sample["delta"]
#             nt = time_sample["nt"]
#             nl = time_sample["nl"]  # npts of left padding
#             nr = time_sample["nr"]  # npts of right padding
#             sem_nt = nt - nl - nr  # number of time sample number in SEM simulation
#             t = np.arange(nt) * dt + (starttime - t0)  # referred to t0
#
#             # ------ get file paths of syn seismograms
#             syn_files = [
#                 "{:s}/{:s}.{:2s}{:1s}{:s}".format(
#                     syn_dir, station_id, syn_band_code, x, syn_suffix
#                 )
#                 for x in syn_orientation_codes
#             ]
#
#             # ------ read in obs, syn seismograms
#             syn_st = read(syn_files[0])
#             syn_st += read(syn_files[1])
#             syn_st += read(syn_files[2])
#
#             # ------ check the same time samples as original syn
#             if not is_equal(
#                 [(tr.stats.starttime, tr.stats.delta, tr.stats.npts) for tr in syn_st]
#             ):
#                 raise Exception(
#                     "%s: not equal time samples in"
#                     " synthetic seismograms." % (station_id)
#                 )
#             tr = syn_st[0]
#
#             if tr.stats.delta != dt:
#                 raise Exception("%s: not the same dt for diff-srcloc!" % (station_id))
#
#             tr_starttime = tr.stats.starttime - nl * dt
#             if tr_starttime != starttime:
#                 raise Exception(
#                     "%s: not the same starttime for diff-srcloc!" % (station_id)
#                 )
#
#             if tr.stats.npts != sem_nt:
#                 raise Exception("%s: not the same npts for diff-srcloc!" % (station_id))
#
#             # ------ read syn seismograms from perturbed source location
#             syn_ENZ = np.zeros((3, nt))
#             for i in range(3):
#                 syn_ENZ[i, nl : (nl + sem_nt)] = syn_st[i].data
#
#             # waveform partial derivatives
#             if "syn" not in waveform:
#                 raise Exception(
#                     "%s: initial syn is not stored in waveform!" % (station_id)
#                 )
#             u0 = waveform["syn"]
#             du = syn_ENZ - u0
#
#             # ------ record derivatives
#             if "waveform_der" not in station:
#                 station["waveform_der"] = {}
#             station["waveform_der"][model_name] = {"du": du}
#
#             # DEBUG: check du
#             # print(dchi
#             # for i in range(3):
#             #  plt.subplot(311+i)
#             #  #plt.plot(t,grn0[i,:],'k', t,syn_ENZ[i,:],'r', t,dg[i,:], 'b')
#             #  plt.plot(t, du[i,:], 'k')
#             # plt.show()
#             if sac_dir:
#                 for i in range(3):
#                     tr.data = du[i, :]
#                     tr.stats.starttime = starttime
#                     tr.stats.delta = dt
#                     tr.stats.npts = nt
#                     out_file = "{:s}/{:s}.{:2s}{:1s}".format(
#                         sac_dir, station_id, syn_band_code, syn_orientation_codes[i]
#                     )
#                     tr.write(out_file, "sac")
#
#     #
#     # ======================================================
#     #
#
#     def measure_hessian_src(
#         self,
#         src_param={"t0": 0, "tau": 0, "xs": 0, "mt": 0},
#         update=True,
#     ):
#         """
#         Approximate Hessian matrix of source parameters based
#         partial derivative waveforms
#
#         Note
#         ----
#         """
#         event = self.data["event"]
#         # src_param = ('dt0','dtau','dxs','dmt')
#         n_srcparam = len(src_param)
#
#         # ------ loop each station
#         station_dict = self.data["station"]
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected statations
#             if station["stat"]["code"] < 1:
#                 continue
#
#             # waveform
#             waveform = station["waveform"]
#             time_sample = waveform["time_sample"]
#             syn_starttime = time_sample["starttime"]
#             syn_delta = time_sample["delta"]
#             syn_nyq = 0.5 / syn_delta
#             syn_nt = time_sample["nt"]
#             syn_nl = time_sample["nl"]
#             syn_nr = time_sample["nr"]
#             syn_times = syn_delta * np.arange(syn_nt)
#             # seismograms
#             obs = waveform["obs"]
#             syn = waveform["syn"]
#
#             # source spectrum (moment-rate function)
#             syn_freq = np.fft.rfftfreq(syn_nt, d=syn_delta)
#             F_src = stf_gauss_spectrum(syn_freq, event["tau"])
#             # spectrum of source derivatives
#             if "t0" in src_param or "tau" in src_param:
#                 F_ds_dt0, F_ds_dtau = stf_gauss_spectrum_der(syn_freq, event["tau"])
#                 # convolve Ds(t)/Dt0,tau with Green's function
#                 du_dt0 = np.fft.irfft(F_ds_dt0 * np.fft.rfft(grn), nt)
#                 du_dtau = np.fft.irfft(F_ds_dtau * np.fft.rfft(grn), nt)
#                 # zero records before origin time (wrap around from the end)
#                 idx = t < -5.0 * tau
#                 du_dt0[:, idx] = 0.0
#                 du_dtau[:, idx] = 0.0
#
#             # waveform derivatives
#             waveform_der = station["waveform_der"]
#
#             # ------ loop each window
#             window_dict = station["window"]
#             for window_id in window_dict:
#                 # window parameters
#                 window = window_dict[window_id]
#                 # skip bad windows
#                 if window["stat"]["code"] < 1:
#                     warnings.warn("Window %s not measured for adj, SKIP" % window_id)
#                     continue
#
#                 # ------ window parameters
#                 # filter
#                 filter_dict = window["filter"]
#                 filter_a = filter_dict["a"]
#                 filter_b = filter_dict["b"]
#                 # taper
#                 win_func = window["taper"]["win"]
#                 # polarity projection
#                 proj_matrix = window["polarity"]["proj_matrix"]
#
#                 # ------ filter obs, syn
#                 # F * d
#                 obs_filt = signal.filtfilt(filter_b, filter_a, obs)
#                 # F * u (u = S*grn)
#                 syn_filt = signal.filtfilt(filter_b, filter_a, grn)
#                 syn_filt = np.fft.irfft(F_src * np.fft.rfft(syn_filt), syn_nt)
#                 # apply window taper and polarity projection
#                 # obs = w * F * d
#                 wFd = np.dot(proj_matrix, obs_filt) * win_func
#                 # syn = w * F * u (u = S*grn)
#                 wFu = np.dot(proj_matrix, syn_filt) * win_func
#                 # norm
#                 norm_wFd = np.sqrt(np.sum(wFd**2))
#                 norm_wFu = np.sqrt(np.sum(wFu**2))
#                 # window normalization factor
#                 Nw = norm_wFd * norm_wFu
#                 # window amplitude raito
#                 Aw = np.sum(wFd * wFu) / norm_wFu**2
#
#                 # DEBUG
#                 # print("Nw: %e %e" % (Nw, window['cc']['Nw'])
#                 # print("Aw: %e %e" % (Aw, window['cc']['Aw'])
#
#                 # ------ filter differential seismograms (w * F * du_dm)
#                 wFdu = {}
#                 for par in src_param:
#                     if par == "t0":
#                         Fdu = signal.filtfilt(filter_b, filter_a, du_dt0)
#                         wFdu[par] = np.dot(proj_matrix, Fdu) * win_func
#                     elif par == "tau":
#                         Fdu = signal.filtfilt(filter_b, filter_a, du_dtau)
#                         wFdu[par] = np.dot(proj_matrix, Fdu) * win_func
#                     else:
#                         du = waveform_der[param]["du"]
#                         Fdu = signal.filtfilt(filter_b, filter_a, du)
#                         wFdu[param] = np.dot(proj_matrix, Fdu) * win_func
#
#                 # ------ hessian src
#                 # chi: zero-lag correlation coef. between wFu and wFd
#                 # hessian: ddchi_dmdm
#                 hessian_src = {}
#                 for i in range(n_srcparam):
#                     for j in range(i, n_srcparam):
#                         par1 = src_param[i]
#                         par2 = src_param[j]
#                         wFdu1 = wFdu[par1]
#                         wFdu2 = wFdu[par2]
#                         wFdu1_wFdu2 = np.sum(wFdu1 * wFdu2)
#                         wFu_wFdu1 = np.sum(wFu * wFdu1)
#                         wFu_wFdu2 = np.sum(wFu * wFdu2)
#                         wFd_wFdu1 = np.sum(wFd * wFdu1)
#                         wFd_wFdu2 = np.sum(wFd * wFdu2)
#                         key12 = (par1, par2)
#                         hessian_src[key12] = (
#                             -Aw * wFdu1_wFdu2
#                             + (
#                                 3.0 * Aw * wFu_wFdu1 * wFu_wFdu2
#                                 - wFu_wFdu1 * wFd_wFdu2
#                                 - wFu_wFdu2 * wFd_wFdu1
#                             )
#                             / norm_wFu**2
#                         ) / Nw
#
#                 # ------ record results
#                 if "hessian_src" not in window:
#                     window["hessian_src"] = hessian_src
#                     window["stat"] = {
#                         "code": 2,
#                         "msg": "add hessian_src on " + UTCDateTime.now().isoformat(),
#                     }
#                 elif update:
#                     window["hessian_src"].update(hessian_src)
#                     window["stat"] = {
#                         "code": 2,
#                         "msg": "update hessian_src on " + UTCDateTime.now().isoformat(),
#                     }
#                 else:
#                     warnings.warn("hessian_src already set, nothing changed")
#             # end for window_id in windows:
#         # endfor station_id in station_dict:
#
#     # enddef measure_windows_for_one_station(self,
#
#     #
#     # ======================================================
#     #
#
#     # def update_source(self, src_param=):
#     #   """ Update source parameters based on waveform derivatives and hessian
#     #   """
#     #   event = self.data['event']
#     #   src_param = ('dt0','dtau','dxs','dmt')
#     #   n_srcparam = len(src_param)
#
#     #   dchi_dm = np.zeros(n_srcparam)
#     #   hessian = np.zeros([n_srcparam,n_srcparam])
#
#     #   #------ get dchi_dm and Hessian
#     #   #-- loop each station
#     #   station_dict = self.data['station']
#     #   for station_id in station_dict:
#     #     station = station_dict[station_id]
#     #     # skip rejected statations
#     #     if station['stat']['code'] < 0:
#     #       continue
#
#     #     # dchi_dm
#     #     for i in range(n_srcparam):
#     #       key = src_param[i]
#     #       dchi_dm[i] += station['waveform_der'][key]['dchi']
#
#     #     #-- loop each window
#     #     window_dict = station['window']
#     #     for window_id in window_dict:
#     #       # window parameters
#     #       window = window_dict[window_id]
#     #       # skip bad windows
#     #       if window['stat']['code'] < 1:
#     #         warnings.warn("Window %s not measured for adj, SKIP" % window_id)
#     #         continue
#     #       #
#     #       weight = window['weight']
#     #       hessian_win = window['hessian_src']
#     #       for i in range(n_srcparam):
#     #         for j in range(i, n_srcparam):
#     #           par1 = src_param[i]
#     #           par2 = src_param[j]
#     #           key = (par1,par2)
#     #           hessian[i,j] += weight * hessian_win[key]
#     #     #end for window_id in windows:
#
#     #   #end for station_id in station_dict:
#
#     #   for i in range(n_srcparam):
#     #     for j in range(i+1, n_srcparam):
#     #         hessian[j,i] = hessian[i,j]
#
#     #   print("dchi_dm:"
#     #   print(dchi_dm
#
#     #   print("hessian:"
#     #   print(hessian
#
#     #   print("====== 0:4:"
#     #   w, v = np.linalg.eigh(hessian, UPLO='U')
#     #   print(w
#     #   print(v
#     #   x, residual, rank, sigval = np.linalg.lstsq(hessian, -dchi_dm)
#     #   print(" inv(hessian)*(-1.0 * dchi_dm): \n", x
#     #   print("dt0: \n", x[0]
#     #   print("dtau:\n", x[1]
#     #   print("dxs: \n", x[2]*self.data['src_perturb']['xs']
#     #   print("dmt: \n", x[3]*self.data['src_perturb']['mt']
#
#     #   print("====== only 0:3"
#     #   h3 = hessian[0:3,0:3]
#     #   v3 = dchi_dm[0:3]
#     #   w, v = np.linalg.eigh(h3, UPLO='U')
#     #   print(w
#     #   print(v
#     #   x, residual, rank, sigval = np.linalg.lstsq(h3, -v3)
#     #   print("inv(hessian)*(-1.0 * dchi_dm): \n", x
#     #   print("dt0: \n", x[0]
#     #   print("dtau:\n", x[1]
#     #   print("dxs: \n", x[2]*self.data['src_perturb']['xs']
#     #   #print("dmt: \n", x[3]*self.data['src_perturb']['dmt']
#
#     #   print("====== only 0:2"
#     #   h3 = hessian[0:2,0:2]
#     #   v3 = dchi_dm[0:2]
#     #   w, v = np.linalg.eigh(h3, UPLO='U')
#     #   print(w
#     #   print(v
#     #   x, residual, rank, sigval = np.linalg.lstsq(h3, -v3)
#     #   print("inv(hessian)*(-1.0 * dchi_dm): \n", x
#     #   print("dt0: \n", x[0]
#     #   print("dtau:\n", x[1]
#
#     #   print("====== only 0,2"
#     #   idx = [0,2]
#     #   hess = hessian[idx,:][:,idx]
#     #   kernel = dchi_dm[idx]
#     #   print(hess
#     #   eigs, eigv = np.linalg.eigh(hess, UPLO='U')
#     #   print(eigs
#     #   print(eigv
#     #   x, residual, rank, sigval = np.linalg.lstsq(hess, -kernel)
#     #   print("inv(hessian)*(-1.0 * dchi_dm): \n", x
#     #   print("dt0: \n", x[0]
#     #   print("dxs:\n", x[1]
#
#     # #enddef measure_windows_for_one_station(self,
#
#     #
#     # ======================================================
#     #
#
#     def cc_linearized_seismogram_for_source(
#         self,
#         dm={"t0": None, "xs_mt": None},
#         plot=False,
#     ):
#         """
#         calculate normalized zero-lag cc for perturbed seismograms
#         by linear combination of waveform derivatives of source parameters (t0,tau,xs,M)
#
#         dm={<model_name>:<model_step_sizes>, ...}
#
#         return wcc_sum, weight_sum
#         """
#         # check model vectors in dm have the same length
#         if (not dm) or len(dm) == 0:
#             error_str = "dm must not be empty"
#             raise Exception(error_str)
#         model_num = len(dm)
#
#         nstep = []
#         for model_name in dm:
#             if dm[model_name].ndim != 1:
#                 error_str = "dm[%s] must be vector" % model_name
#                 raise Exception(error_str)
#             nstep.append(np.size(dm[model_name]))
#         if not is_equal(nstep) or nstep[0] < 1:
#             error_str = "vectors in dm must have the same non-zero length"
#             raise Exception(error_str)
#         nstep = nstep[0]
#
#         # check parameters
#         event = self.data["event"]
#         if "tau" in dm:
#             tau = dm["tau"] + event["tau"]
#             if any(tau < 0):
#                 error_str = (
#                     "dm['dtau'] has invalid values (event['tau']=%f)!" % event["tau"]
#                 )
#                 raise Exception(error_str)
#
#         # ------ loop each station
#         # sum of weighted normalized zero-lag cc at each model grid
#         wcc_sum = np.zeros(nstep)
#         # sum of all windows' weighting
#         weight_sum = 0.0
#
#         station_dict = self.data["station"]
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected stations
#             if station["stat"]["code"] < 0:
#                 continue
#
#             # print(station_id, flush=True)
#             # print(datetime.datetime.utcnow(), flush=True)
#
#             # check if model parameter included in waveform_der
#             waveform_der = station["waveform_der"]
#             for model_name in dm:
#                 if (model_name not in ["t0", "tau"]) and (
#                     model_name not in waveform_der
#                 ):
#                     error_str = "%s not in waveform_der of %s" % (
#                         model_name,
#                         station_id,
#                     )
#                     raise Exception(error_str)
#
#             # ---- get seismograms: obs,grn
#             waveform = station["waveform"]
#             obs = waveform["obs"]
#             grn = waveform["grn"]
#             # time samples
#             time_sample = waveform["time_sample"]
#             syn_starttime = time_sample["starttime"]
#             syn_delta = time_sample["delta"]
#             syn_nt = time_sample["nt"]
#             syn_nl = time_sample["nl"]
#             syn_nr = time_sample["nr"]
#             syn_freq = np.fft.rfftfreq(syn_nt, d=syn_delta)
#
#             # ---- measure misfit
#             window_dict = station["window"]
#             for window_id in window_dict:
#                 window = window_dict[window_id]
#                 # skip bad windows
#                 if window["stat"]["code"] < 1:
#                     warnings.warn("Window %s not measured for adj, SKIP" % window_id)
#                     continue
#                 # window weight
#                 weight = window["weight"]
#                 # skip window with zero weight
#                 if np.isclose(weight, 0.0):
#                     continue
#                 weight_sum += weight
#                 # filter
#                 filter_dict = window["filter"]
#                 filter_a = filter_dict["a"]
#                 filter_b = filter_dict["b"]
#                 # taper
#                 win_func = window["taper"]["win"]
#                 win_starttime = window["taper"]["starttime"]
#                 win_endtime = window["taper"]["endtime"]
#                 if plot:
#                     syn_times = np.arange(syn_nt) * syn_delta
#                     plt.plot(syn_times, win_func)
#                     plt.show()
#                     plt.title("wind_func")
#                 # polarity projection
#                 proj_matrix = window["polarity"]["proj_matrix"]
#                 # -- filter,project,taper obs
#                 # F * d
#                 obs_filt = signal.filtfilt(filter_b, filter_a, obs)
#                 # w * p * F * d (window,project,filter)
#                 wpFd = np.dot(proj_matrix, obs_filt) * win_func
#                 norm_wpFd = np.sqrt(np.sum(wpFd**2))
#                 # -- filter,project grn
#                 # F * g
#                 grn_filt = signal.filtfilt(filter_b, filter_a, grn)
#                 # p * F * g
#                 pFg = np.dot(proj_matrix, grn_filt)
#                 # if plot:
#                 #  F_src = stf_gauss_spectrum(syn_freq, event['tau'])
#                 #  # S * F * g
#                 #  syn_filt = np.fft.irfft(F_src*np.fft.rfft(grn_filt), syn_nt)
#                 #  # w * p * S * F * g
#                 #  wpFu = np.dot(proj_matrix, syn_filt) * win_func
#                 # -- filter,project dg: pFdg
#                 pFdg = {}
#                 for model_name in dm:
#                     # exclude source time function parameters
#                     if model_name not in ["t0", "tau"]:
#                         dg = waveform_der[model_name]["dg"]
#                         dg_filt = signal.filtfilt(filter_b, filter_a, dg)
#                         pFdg[model_name] = np.dot(proj_matrix, dg_filt)
#                 # -- misfit function: zero-lag cc
#                 for istep in range(nstep):
#                     # perturbed grn: pFg1
#                     pFg1 = np.zeros((3, syn_nt))
#                     pFg1 += pFg
#                     for model_name in dm:
#                         # exclude source time function
#                         if model_name not in ["t0", "tau"]:
#                             pFg1 += dm[model_name][istep] * pFdg[model_name]
#                     # perturbed source time function
#                     dt0 = 0.0
#                     if "t0" in dm:
#                         dt0 = dm["t0"][istep]
#                     dtau = 0.0
#                     if "tau" in dm:
#                         dtau = dm["tau"][istep]
#                     F_src = stf_gauss_spectrum(syn_freq, event["tau"] + dtau)
#                     # perturbed syn: w * S * p * F * g1
#                     phase_shift = np.exp(-2.0j * np.pi * syn_freq * dt0)
#                     wpFu1 = (
#                         np.fft.irfft(phase_shift * F_src * np.fft.rfft(pFg1), syn_nt)
#                         * win_func
#                     )
#                     norm_wpFu1 = np.sqrt(np.sum(wpFu1**2))
#                     Nw = norm_wpFd * norm_wpFu1
#                     # normalized cc between obs and perturbed syn
#                     cc_wpFd_wpFu1 = np.sum(wpFd * wpFu1) / Nw
#                     # weighted cc
#                     wcc_sum[istep] += weight * cc_wpFd_wpFu1
#                     # DEBUG
#                     # if plot:
#                     #  syn_npts = syn_nt - syn_nl - syn_nr
#                     #  syn_orientation_codes = ['E', 'N', 'Z']
#                     #  syn_times = np.arange(syn_nt) * syn_delta
#                     #  Amax_obs = np.max(np.sum(wpFd**2, axis=0))**0.5
#                     #  Amax_syn = np.max(np.sum(wpFu**2, axis=0))**0.5
#                     #  Amax_syn1 = np.max(np.sum(wpFu1**2, axis=0))**0.5
#                     #  win_b = win_starttime - syn_starttime
#                     #  win_e = win_endtime - syn_starttime
#                     #  obs_filt_proj = np.dot(proj_matrix, obs_filt)
#                     #  syn_filt_proj = np.dot(proj_matrix, syn_filt)
#                     #  for i in range(3):
#                     #    plt.subplot(311+i)
#                     #    if i == 0:
#                     #      title_str = "%s.%s " % (station_id, window_id)
#                     #      for model_name in dm:
#                     #        title_str += "%s:%.2f " \
#                     #            % (model_name, dm[model_name][idx_model])
#                     #      title_str += "NCCwin:%.2f" % (cc_wpFd_wpFu1)
#                     #      plt.title(title_str)
#
#                     #    idx_plt = range(syn_nl,(syn_nl+syn_npts))
#                     #    # whole trace, projected
#                     #    plt.plot(syn_times[idx_plt], obs_filt_proj[i,idx_plt]/Amax_obs,
#                     #        'k', linewidth=0.2)
#                     #    plt.plot(syn_times[idx_plt], syn_filt_proj[i,idx_plt]/Amax_syn,
#                     #        'b', linewidth=0.2)
#                     #    # windowed trace
#                     #    idx_plt = (win_b <= syn_times) & (syn_times <= win_e)
#                     #    plt.plot(syn_times[idx_plt], wpFd[i,idx_plt]/Amax_obs,
#                     #        'k', linewidth=1.0)
#                     #    plt.plot(syn_times[idx_plt], wpFu[i,idx_plt]/Amax_syn,
#                     #        'b', linewidth=1.0)
#                     #    plt.plot(syn_times[idx_plt], wpFu1[i,idx_plt]/Amax_syn1,
#                     #        'r', linewidth=1.0)
#
#                     #    plt.xlim((syn_times[syn_nl], syn_times[syn_nl+syn_npts-1]))
#                     #    plt.ylim((-1.5, 1.5))
#                     #    plt.ylabel(syn_orientation_codes[i])
#                     #  plt.show()
#
#             # end for window_id in window_dict:
#         # end for station_id in station_dict:
#
#         return wcc_sum, weight_sum
#
#     #
#     # ======================================================
#     #
#
#     #  def search1d_cc_perturbed_seismogram(self,
#     #      dm_range = {
#     #        't0': [-5,5],
#     #        'tau': [-5,5],
#     #        'xs': [-5,5],
#     #        'mt': [-5,5],
#     #        },
#     ##     dist_min=0.0,
#     ##     dist_max=180.0,
#     #      ngrid=4,
#     #      max_niter=5,
#     #      range_ratio=0.7,
#     #      plot_seism=False,
#     #      log_file="search1d_cc.log",
#     #      cmt_file="CMTSOLUTION.search1d",
#     #      ):
#     #    """
#     #    calculate misfit over 2D model grids based on perturbed seismograms
#     #
#     #    Parameters
#     #    ----------
#     #    range_ratio : scalar
#     #      between (0, 1), used to shrink the search range to cc values above range_ratio
#     #      between cc_max and cc_min.
#     #
#     #    Return
#     #    ------
#     #    optimal dm
#     #
#     #    """
#     #    # check validity of dm_range
#     #    event = self.data['event']
#     #    tau = event['tau']
#     #    if 'tau' in dm_range:
#     #      tau_lim = dm_range['tau']
#     #      if (min(tau_lim) + tau) < 1.0:
#     #        warnings.warn("Too small dtau value, increase to avoid negative tau")
#     #        dm_range['tau'][0] = -0.95*event['tau']
#     #        print(dm_range)
#     #
#     #    # grid parameters
#     #    model_name = [x for x in dm_range]
#     #
#     #    range_ratio = float(range_ratio)
#     #    if range_ratio < 0.5 or range_ratio > 0.9:
#     #      raise ValueError("range_ratio should be within (0.5, 0.9)")
#     #
#     #    dm_opt = {}
#     #    for par in model_name:
#     #      dm_opt[par] = 0.0
#     #
#     #    # prepare log file
#     #    fid = open(log_file,"w")
#     #
#     #    niter = 0
#     #    while niter < max_niter:
#     #      fid.write("====== iteration: %d\n" % (niter))
#     #      fid.flush()
#     #
#     #      # loop each parameter
#     #      for par in model_name:
#     #        fid.write("------ %s\n" % (par))
#     #        fid.flush()
#     #        xlim = np.array(dm_range[par])
#     #
#     #        if xlim.size == 1:
#     #          dm_opt[par] = xlim[0]
#     #          fid.write("fixed at %.5e\n" % (dm_opt[par]))
#     #          fid.flush()
#     #          continue
#     #        else:
#     #          fid.write("range: %.5e, %.5e\n" % (np.min(xlim), np.max(xlim)))
#     #          fid.flush()
#     #          xx = np.linspace(np.min(xlim), np.max(xlim), ngrid)
#     #
#     #          # search
#     #          dm_grid = {}
#     #          dm_grid[par] = xx
#     #          for par_fix in model_name:
#     #            if par_fix != par:
#     #              dm_grid[par_fix] = np.ones(ngrid) * dm_opt[par_fix]
#     #          cc, weight = self.cc_perturbed_seismogram(
#     #              dm=dm_grid,
#     ##             dist_min=dist_min,
#     ##             dist_max=dist_max,
#     #              plot=plot_seism)
#     #          cc /= weight # weighted average of normalized zero-lag CC
#     #
#     #          # interpolation
#     #          x_i = np.linspace(np.min(xlim), np.max(xlim), 100)
#     #          interp = interpolate.interp1d(xx, cc, kind='cubic')
#     #          cc_i = interp(x_i)
#     #          imax = np.argmax(cc_i)
#     #          x_max = x_i[imax]
#     #          cc_max = cc_i[imax]
#     #
#     #          #change search range
#     #          imin = np.argmin(cc_i)
#     #          cc_min = cc_i[imin]
#     #          cc_thred = cc_min + (cc_max-cc_min)*range_ratio
#     #          x1 = x_i[cc_i>cc_thred]
#     #          dm_range[par] = [np.min(x1), np.max(x1)]
#     #
#     #          dm_opt[par] = x_max
#     #          fid.write("maximum at %.5e with cc %.5e\n" % (dm_opt[par], cc_max))
#     #          fid.flush()
#     #
#     #      # print(out optimal model
#     #      fid.write("------ optimal dm_alpha\n")
#     #      fid.flush()
#     #      for par in model_name:
#     #        fid.write("%s: %.5e\n" %(par, dm_opt[par]))
#     #        fid.flush()
#     #
#     #      niter += 1
#     #    #END while niter < 5:
#     #    fid.close()
#     #
#     #    # make new model
#     #    model = {}
#     #    event = self.data['event']
#     #    src_perturb = self.data['src_perturb']
#     #    for par in dm_opt:
#     #      m0 = event[par]
#     #      dm = src_perturb[par]
#     #      alpha = dm_opt[par]
#     #      model[par] = {
#     #          'm0':m0,
#     #          'dm':dm,
#     #          'alpha':alpha,
#     #          'm1':m0+alpha*dm
#     #          }
#     #
#     #    # write out new CMTSOLUTION file
#     #    t0 = event['t0']
#     #    if 't0' in model: t0 = model['t0']['m1']
#     #    # modify origin time in header line to have centroid time
#     #    header = event['header'].split()
#     #    header[1] = "{:04d}".format(t0.year)
#     #    header[2] = "{:02d}".format(t0.month)
#     #    header[3] = "{:02d}".format(t0.day)
#     #    header[4] = "{:02d}".format(t0.hour)
#     #    header[5] = "{:02d}".format(t0.minute)
#     #    header[6] = "{:07.4f}".format(t0.second + 1.0e-6*t0.microsecond)
#     #
#     #    tau = event['tau']
#     #    if 'tau' in model: tau = model['tau']['m1']
#     #    xs = event['xs']
#     #    if 'xs' in model: xs = model['xs']['m1']
#     #    mt = event['mt']
#     #    m0 = (0.5*np.sum(mt**2))**0.5
#     #    if 'mt' in model: mt = model['mt']['m1']
#     #    # zero trace
#     #    mt = mt - np.identity(3)*np.trace(mt)/3.0
#     #    # scale to orignal moment
#     #    mt = mt/(0.5*np.sum(mt**2))**0.5
#     #    mt *= m0
#     #
#     #    with open(cmt_file, 'w') as fp:
#     #      fp.write('%s \n' % (' '.join(header)))
#     #      fp.write('%-18s %s\n' % ('event name:', event['id']))
#     #      fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
#     #      fp.write('%-18s %+15.8E\n' % ('tau(s):',   tau))
#     #      fp.write('%-18s %+15.8E\n' % ('x(m):',     xs[0]))
#     #      fp.write('%-18s %+15.8E\n' % ('y(m):',     xs[1]))
#     #      fp.write('%-18s %+15.8E\n' % ('z(m):',     xs[2]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mt[0,0]))
#     #      fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', mt[1,1]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mt[2,2]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mt[0,1]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mt[0,2]))
#     #      fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', mt[1,2]))
#
#     #
#     # ======================================================
#     #
#
#     #  def grid_cc_perturbed_seismogram(self,
#     #      dm = {
#     #        't0': np.linspace(-10,0,11),
#     #        'xs': np.linspace(-5,5,11),
#     #        },
#     #      axes=[ ('t0','xs'), ('t0',), ('xs',) ],
#     #      outfig="grid_cc.pdf",
#     #      plot_seism=False,
#     #      cmt_file="CMTSOLUTION.grid_cc"
#     #      ):
#     #    """ calculate misfit over 2D model grids based on perturbed seismograms
#     #    """
#     #    # grid parameters
#     #    model_num = len(dm)
#     #    model_name = [x for x in dm]
#     #
#     #    # check parameters
#     #    for xy in axes:
#     #      if len(xy) < 1 or len(xy) > 2:
#     #        error_str = "axes should have one or two elements"
#     #        raise Exception(error_str)
#     #      for axis in xy:
#     #        if axis not in model_name:
#     #          error_str = "axis(%s) not in dm" % axis
#     #          raise Exception(error_str)
#     #
#     #    # model grid
#     #    x = [np.array(dm[s]) for s in model_name]
#     #    xx = np.meshgrid(*x, indexing='ij')
#     #    nn = xx[0].size
#     #
#     #    # calculate cc values over all grid points
#     #    dm_grid = {}
#     #    for i in range(model_num):
#     #      par = model_name[i]
#     #      dm_grid[par] = xx[i].flatten()
#     #
#     #    zz, weight = self.cc_perturbed_seismogram(dm=dm_grid, plot=plot_seism)
#     #    zz /= weight # weighted average of normalized zero-lag CC
#     #    zz = zz.reshape(xx[0].shape)
#     #
#     #    # get maximum cc value
#     #    imax = np.argmax(zz)
#     #    idx_max = np.unravel_index(imax, xx[0].shape)
#     #    zz_max = zz[idx_max]
#     #
#     #    zz_min = np.min(zz)
#     #
#     #    # plot out results
#     #    with PdfPages(outfig) as pdf:
#     #      for xy in axes:
#     #        # plot 1D section
#     #        if (len(xy) == 1) or (xy[0] == xy[1]):
#     #          ix = model_name.index(xy[0])
#     #          idx = [range(len(v)) for v in x]
#     #          for i in range(model_num):
#     #            if i != ix:
#     #              idx[i] = [idx_max[i],]
#     #          zz_1d = zz[np.ix_(*idx)].squeeze()
#     #
#     #          xx_1d = x[ix]
#     #          x_i = np.linspace(np.min(xx_1d), np.max(xx_1d), 100)
#     #
#     #          interp = interpolate.interp1d(xx_1d, zz_1d, kind='cubic')
#     #          zz_i = interp(x_i)
#     #          imax_i = np.argmax(zz_i)
#     #          x_i_max = x_i[imax_i]
#     #          zz_i_max = zz_i[imax_i]
#     #
#     #          # plot cross-sections through maximum point
#     #          #fig = plt.figure(figsize=(8.5,11))
#     #          fig = plt.figure()
#     #          plt.plot(xx_1d, zz_1d, 'ro', x_i, zz_i)
#     #
#     #          title_str = "Weighted average CC  (%.2e, %.2e)" % (x_i_max, zz_i_max)
#     #          plt.title(title_str)
#     #
#     #          xlabel_str = "alpha (* d%s)" % model_name[ix]
#     #          plt.xlabel(xlabel_str)
#     #          plt.ylabel("weighted avg. CC")
#     #
#     #          pdf.savefig()
#     #          plt.close()
#     #
#     #        # plot 2D cross section
#     #        else:
#     #          ix = model_name.index(xy[0])
#     #          iy = model_name.index(xy[1])
#     #          idx = [range(len(v)) for v in x]
#     #          for i in range(model_num):
#     #            if i != ix and i != iy:
#     #              idx[i] = [idx_max[i],]
#     #          idx = np.ix_(*idx)
#     #          zz_2d = zz[idx].squeeze()
#     #          xx_2d = xx[ix][idx].squeeze()
#     #          yy_2d = xx[iy][idx].squeeze()
#     #          print(xx_2d)
#     #          print(yy_2d)
#     #          print(zz_2d)
#     #          if ix > iy:
#     #            zz_2d = np.transpose(zz_2d)
#     #            xx_2d = np.transpose(xx_2d)
#     #            yy_2d = np.transpose(yy_2d)
#     #
#     #          # make interpolation grid
#     #          x_i = np.linspace(np.min(xx_2d), np.max(xx_2d), 100)
#     #          y_i = np.linspace(np.min(yy_2d), np.max(yy_2d), 100)
#     #          xx_i, yy_i = np.meshgrid(x_i, y_i)
#     #
#     #          interp = interpolate.interp2d(xx_2d, yy_2d, zz_2d, kind='cubic')
#     #          zz_i = interp(x_i, y_i)
#     #          imax_i = np.argmax(zz_i)
#     #          idx_max_i = np.unravel_index(imax_i, xx_i.shape)
#     #          zz_i_max = zz_i[idx_max_i]
#     #          x_i_max = xx_i[idx_max_i]
#     #          y_i_max = yy_i[idx_max_i]
#     #
#     #          # plot 2D surface through the maximum
#     #          #fig = plt.figure(figsize=(8.5, 11))
#     #          fig = plt.figure()
#     #
#     #          levels = zz_min + (zz_max-zz_min)*np.linspace(0.5, 1.0, 20)
#     #          #cs = plt.contour(xx_2d, yy_2d, zz_2d, levels, colors='k')
#     #          cs = plt.contour(xx_i, yy_i, zz_i, levels, colors='k')
#     #          plt.clabel(cs, fontsize=9, inline=1)
#     #
#     #          #x_max = xx[ix][idx_max]
#     #          #y_max = xx[iy][idx_max]
#     #          #text_str = "(%.1f,%.1f)" % (x_max, y_max)
#     #          #plt.text(x_max, y_max, text_str,
#     #          #    horizontalalignment="center", verticalalignment="top")
#     #          #plt.plot(x_max, y_max, 'ko')
#     #
#     #          #text_str = "%.2e,%.2e" % (x_i_max, y_i_max)
#     #          #plt.text(x_i_max, y_i_max, text_str,
#     #          #    horizontalalignment="center", verticalalignment="top")
#     #          plt.plot(x_i_max, y_i_max, 'k+')
#     #
#     #          title_str = "Weighted average CC (%.2e,%.2e)" % (x_i_max, y_i_max)
#     #          plt.title(title_str)
#     #
#     #          xlabel_str = "alpha (* d%s)" % model_name[ix]
#     #          plt.xlabel(xlabel_str)
#     #          ylabel_str = "alpha (* d%s)" % model_name[iy]
#     #          plt.ylabel(ylabel_str)
#     #
#     #          pdf.savefig()
#     #          plt.close()
#     #
#     #    # get best model
#     #    event = self.data['event']
#     #    src_perturb = self.data['src_perturb']
#     #
#     #    model = {}
#     #    for par in dm:
#     #      m0 = event[par]
#     #      dm = src_perturb[par]
#     #      ix = model_name.index(par)
#     #      alpha = xx[ix][idx_max]
#     #      model[par] = {
#     #          'm0':m0,
#     #          'dm':dm,
#     #          'alpha':alpha,
#     #          'm1':m0+alpha*dm
#     #          }
#     #
#     #    #------ write out new CMTSOLUTION file
#     #    t0 = event['t0']
#     #    if 't0' in model: t0 = model['t0']['m1']
#     #    # modify origin time in header line to have centroid time
#     #    header = event['header'].split()
#     #    header[1] = str(t0.year)
#     #    header[2] = str(t0.month)
#     #    header[3] = str(t0.day)
#     #    header[4] = str(t0.hour)
#     #    header[5] = str(t0.minute)
#     #    header[6] = str(t0.second + 1.0e-6*t0.microsecond)
#     #
#     #    tau = event['tau']
#     #    if 'tau' in model: tau = model['tau']['m1']
#     #    xs = event['xs']
#     #    if 'xs' in model: xs = model['xs']['m1']
#     #    mt = event['mt']
#     #    if 'mt' in model: mt = model['mt']['m1']
#     #
#     #    with open(cmt_file, 'w') as fp:
#     #      fp.write('%s | wCCmax %f weight %f\n'
#     #          % (' '.join(header), zz_max, weight))
#     #      fp.write('%-18s %s\n' % ('event name:', "grid_cc"))
#     #      fp.write('%-18s %+15.8E\n' % ('t0(s):',    0.0))
#     #      fp.write('%-18s %+15.8E\n' % ('tau(s):',   tau))
#     #      fp.write('%-18s %+15.8E\n' % ('x(m):',     xs[0]))
#     #      fp.write('%-18s %+15.8E\n' % ('y(m):',     xs[1]))
#     #      fp.write('%-18s %+15.8E\n' % ('z(m):',     xs[2]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mxx(N*m):', mt[0,0]))
#     #      fp.write('%-18s %+15.8E\n' % ('Myy(N*m):', mt[1,1]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mzz(N*m):', mt[2,2]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mxy(N*m):', mt[0,1]))
#     #      fp.write('%-18s %+15.8E\n' % ('Mxz(N*m):', mt[0,2]))
#     #      fp.write('%-18s %+15.8E\n' % ('Myz(N*m):', mt[1,2]))
#     #
#     #    return model, zz_max, weight
#
#     #
#     # ======================================================
#     #
#
#     #  def relocate_1d(self,
#     #      event_id,
#     #      window_id_list=['F.p,P', 'F.s,S'],
#     #      fix_depth=False,
#     #      out_cmt_file=None):
#     #    """relocate event using ray path in reference earth model
#     #    """
#     #    # check inputs
#     #    events = self.data['events']
#     #    if event_id not in events:
#     #      print("[ERROR] %s does NOT exist. Exit" % (event_id))
#     #      sys.exit()
#     #
#     #    # select windows
#     #    sta_win_id_list = []
#     #    event = events[event_id]
#     #    stations = event['stations']
#     #    for station_id in stations:
#     #
#     #      station = stations[station_id]
#     #      if station['stat']['code'] < 0:
#     #        continue
#     #
#     #      windows = station['windows']
#     #      for window_id in window_id_list:
#     #        if window_id not in windows:
#     #          continue
#     #
#     #        window = windows[window_id]
#     #        if window['stat']['code'] != 1:
#     #          continue
#     #
#     #        misfit = window['misfit']
#     #        #if window['quality']['SNR'] < min_SNR or \
#     #        #    misfit['CC0'] < min_CC0 or \
#     #        #    misfit['CCmax'] < min_CCmax or\
#     #        #    abs(misfit['CC_time_shift']) > max_CC_time_shift:
#     #        #  continue
#     #
#     #        sta_win_id = (station_id, window_id)
#     #        sta_win_id_list.append(sta_win_id)
#     #
#     #    # create sensitivity matrix G in local NED coordinate
#     #    # G * dm  = dt_cc
#     #    # G: [[-px_1, -py_1, -pz_1, 1.0], # ray1
#     #    #   [-px_2, -py_2, -pz_2, 1.0], # ray2
#     #    #   ...]
#     #    # dm: [dNorth(km), dEast, dDepth, dT(sec)]
#     #    # dt_cc: [dt1, dt2, ...]
#     #    n = len(sta_win_id_list)
#     #    G = np.zeros((n, 4))
#     #    dt_cc = np.zeros(n)
#     #    R_Earth_km = 6371.0
#     #    gcmt = event['gcmt']
#     #    evdp = gcmt['depth']
#     #    for i in range(n):
#     #      sta_win_id = sta_win_id_list[i]
#     #      station_id = sta_win_id[0]
#     #      window_id = sta_win_id[1]
#     #
#     #      station = stations[station_id]
#     #      meta = station['meta']
#     #      window = station['windows'][window_id]
#     #      phase = window['phase']
#     #      misfit = window['misfit']
#     #      weight = window['weight']
#     #
#     #      azimuth = np.deg2rad(meta['azimuth'])
#     #      takeoff_angle = phase['takeoff_angle']
#     #      takeoff = np.deg2rad(takeoff_angle + 180.0*(takeoff_angle<0))
#     #      ray_param = phase['ray_param']
#     #      slowness = ray_param / (R_Earth_km - evdp) #unit: s/km
#     #      # local coordinate: NED
#     #      pd = np.cos(takeoff) * slowness
#     #      pn = np.cos(azimuth) * np.sin(takeoff) * slowness
#     #      pe = np.sin(azimuth) * np.sin(takeoff) * slowness
#     #      # create sensitivity matrix
#     #      G[i,:] = weight * np.array([-pn, -pe, -pd, 1.0]) # -p: from receiver to source
#     #      dt_cc[i] = weight * misfit['CC_time_shift']
#     #
#     #    #linearized inversion (can be extended to second order using dynamic ray-tracing)
#     #    if fix_depth:
#     #      G[:, 2] = 0.0
#     #    dm, residual, rank, sigval = np.linalg.lstsq(G, dt_cc)
#     #
#     #    # convert dm from NED to ECEF coordinate
#     #    evla = gcmt['latitude']
#     #    evlo = gcmt['longitude']
#     #
#     #    slat = np.sin(np.deg2rad(evla))
#     #    clat = np.cos(np.deg2rad(evla))
#     #    slon = np.sin(np.deg2rad(evlo))
#     #    clon = np.cos(np.deg2rad(evlo))
#     #
#     #    N = [-slat*clon, -slat*slon, clat]
#     #    E = [-slon, clon, 0.0]
#     #    D = [-clat*clon, -clat*slon, -slat]
#     #
#     #    ev_dx = N[0]*dm[0] + E[0]*dm[1] + D[0]*dm[2]
#     #    ev_dy = N[1]*dm[0] + E[1]*dm[1] + D[1]*dm[2]
#     #    ev_dz = N[2]*dm[0] + E[2]*dm[1] + D[2]*dm[2]
#     #    ev_dt = dm[3]
#     #
#     #    # initialize pyproj objects
#     #    #geod = pyproj.Geod(ellps='WGS84')
#     #    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
#     #    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
#     #
#     #    # old location in ECEF (meters)
#     #    evx, evy, evz = pyproj.transform(lla, ecef, evlo, evla, -1000.0*evdp)
#     #
#     #    # new location in ECEF (meters)
#     #    evx1 = evx + ev_dx*1000.0
#     #    evy1 = evy + ev_dy*1000.0
#     #    evz1 = evz + ev_dz*1000.0
#     #    # in LLA
#     #    evlo1, evla1, evalt1 = pyproj.transform(ecef, lla, evx1, evy1, evz1)
#     #    evdp1 = -evalt1/1000.0
#     #
#     #    # residuals
#     #    # linearized modelling
#     #    dt_syn = G.dot(dm)
#     #    dt_res = dt_cc - dt_syn
#     #
#     #    # make results
#     #    new_centroid_time = UTCDateTime(gcmt['centroid_time']) + ev_dt
#     #    reloc_dict = {
#     #        'window_id_list': window_id_list,
#     #        'singular_value': sigval.tolist(),
#     #        'dm': {'dNorth':dm[0], 'dEast':dm[1], 'dDepth':dm[2],
#     #          'dT':dm[3]},
#     #        'latitude':evla1,
#     #        'longitude':evlo1,
#     #        'depth':evdp1,
#     #        'centroid_time': str(new_centroid_time),
#     #        'data': {'num':n, 'mean':np.mean(dt_cc), 'std':np.std(dt_cc)},
#     #        'residual': {'mean':np.mean(dt_res), 'std':np.std(dt_res)} }
#     #
#     #    event['relocate'] = reloc_dict
#     #
#     #    # make new CMTSOLUTION file
#     #    if out_cmt_file:
#     #      M = gcmt['moment_tensor']
#     #      with open(out_cmt_file, 'w') as fp:
#     #        # header line:
#     #        #PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
#     #        # which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region
#     #        fp.write(new_centroid_time.strftime(
#     #          'RELOC %Y %m %d %H %M %S.%f ') + \
#     #          '%.4f %.4f %.1f 0.0 0.0 END\n' % (evla1,evlo1,evdp1) )
#     #        fp.write('event name:    %s\n'   % (event_id))
#     #        fp.write('time shift:    0.0\n'        )
#     #        fp.write('tau:   %.1f\n'   % (gcmt['tau']))
#     #        #fp.write('half duration:   0.0\n'  % (gcmt['tau']))
#     #        fp.write('latitude:    %.4f\n'   % (evla1)   )
#     #        fp.write('longitude:     %.4f\n'   % (evlo1)   )
#     #        fp.write('depth:       %.4f\n'   % (evdp1)   )
#     #        fp.write('Mrr:       %12.4e\n' % (M[0][0]) )
#     #        fp.write('Mtt:       %12.4e\n' % (M[1][1]) )
#     #        fp.write('Mpp:       %12.4e\n' % (M[2][2]) )
#     #        fp.write('Mrt:       %12.4e\n' % (M[0][1]) )
#     #        fp.write('Mrp:       %12.4e\n' % (M[0][2]) )
#     #        fp.write('Mtp:       %12.4e\n' % (M[1][2]) )
#
#     #
#     # ======================================================
#     #
#
#     def cc_linearized_seismogram_for_dmodel(
#         self, dm={"vsv": None, "vsh": None}, plot=False
#     ):
#         """
#         Calculate normalized zero-lag cc between observed and linearized seismograms with waveform_der_dmodel
#         for different step size.
#
#         return weighted_cc_sum[:], weight_sum
#         """
#         # ------ validate inputs
#         # check step length arrays in dm have the same length
#         if (not dm) or len(dm) == 0:
#             error_str = "dm must not be empty"
#             raise Exception(error_str)
#         model_num = len(dm)
#
#         nstep = []
#         for model_name in dm:
#             if dm[model_name].ndim != 1:
#                 error_str = "dm[%s] must be vector" % model_name
#                 raise Exception(error_str)
#             nstep.append(np.size(dm[model_name]))
#         if not is_equal(nstep) or nstep[0] < 1:
#             error_str = "vectors in dm must have the same non-zero length"
#             raise Exception(error_str)
#         nstep = nstep[0]
#
#         # ------ loop each station
#         event = self.data["event"]
#         # sum of weighted normalized zero-lag cc at each model grid
#         wcc_sum = np.zeros(nstep)
#         # sum of all windows' weighting
#         weight_sum = 0.0
#
#         station_dict = self.data["station"]
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected stations
#             if station["stat"]["code"] < 0:
#                 continue
#
#             # check if model parameter included in waveform_der
#             waveform_der = station["waveform_der"]
#             for model_name in dm:
#                 if model_name not in waveform_der:
#                     error_str = "%s: waveform_der['%s'] does not exist." % (
#                         station_id,
#                         model_name,
#                     )
#                     raise Exception(error_str)
#
#             # ---- get seismograms: obs,grn
#             waveform = station["waveform"]
#             obs = waveform["obs"]
#             syn = waveform["syn"]
#             # time samples
#             time_sample = waveform["time_sample"]
#             syn_starttime = time_sample["starttime"]
#             syn_delta = time_sample["delta"]
#             syn_nt = time_sample["nt"]
#             syn_nl = time_sample["nl"]
#             syn_nr = time_sample["nr"]
#             syn_freq = np.fft.rfftfreq(syn_nt, d=syn_delta)
#
#             # ---- measure misfit
#             window_dict = station["window"]
#             for window_id in window_dict:
#                 window = window_dict[window_id]
#                 # skip bad windows
#                 if window["stat"]["code"] < 1:
#                     warnings.warn("Window %s not measured for adj, SKIP" % window_id)
#                     continue
#                 # window weight
#                 weight = window["weight"]
#                 # skip window with zero weight
#                 if np.isclose(weight, 0.0):
#                     continue
#                 weight_sum += weight
#                 # filter
#                 filter_dict = window["filter"]
#                 filter_a = filter_dict["a"]
#                 filter_b = filter_dict["b"]
#                 # taper
#                 win_func = window["taper"]["win"]
#                 win_starttime = window["taper"]["starttime"]
#                 win_endtime = window["taper"]["endtime"]
#                 if plot:
#                     syn_times = np.arange(syn_nt) * syn_delta
#                     plt.plot(syn_times, win_func)
#                     plt.show()
#                     plt.title("wind_func")
#                 # polarity projection
#                 proj_matrix = window["polarity"]["proj_matrix"]
#                 # -- filter,project,taper obs
#                 # F * d
#                 obs_filt = signal.filtfilt(filter_b, filter_a, obs)
#                 # w * p * F * d (window,project,filter)
#                 wpFd = np.dot(proj_matrix, obs_filt) * win_func
#                 norm_wpFd = np.sqrt(np.sum(wpFd**2))
#                 # -- filter,project,taper syn
#                 # F * u
#                 syn_filt = signal.filtfilt(filter_b, filter_a, syn)
#                 # w * p * F * u
#                 wpFu = np.dot(proj_matrix, syn_filt) * win_func
#                 # -- filter,project,taper du: wpFdu
#                 wpFdu = {}
#                 for model_name in dm:
#                     du = waveform_der[model_name]["du"]
#                     du_filt = signal.filtfilt(filter_b, filter_a, du)
#                     # w * p * F * du
#                     wpFdu[model_name] = np.dot(proj_matrix, du_filt) * win_func
#                 # -- misfit function: zero-lag cc
#                 for istep in range(nstep):
#                     wpFu1 = np.zeros((3, syn_nt))
#                     wpFu1 += wpFu
#                     for model_name in dm:
#                         wpFu1 += dm[model_name][istep] * wpFdu[model_name]
#                     norm_wpFu1 = np.sqrt(np.sum(wpFu1**2))
#                     Nw = norm_wpFd * norm_wpFu1
#                     # normalized cc between obs and perturbed syn
#                     cc_wpFd_wpFu1 = np.sum(wpFd * wpFu1) / Nw
#                     # weighted cc
#                     wcc_sum[istep] += weight * cc_wpFd_wpFu1
#                     # DEBUG
#                     if plot:
#                         syn_npts = syn_nt - syn_nl - syn_nr
#                         syn_orientation_codes = ["E", "N", "Z"]
#                         syn_times = np.arange(syn_nt) * syn_delta
#                         Amax_obs = np.max(np.sum(wpFd**2, axis=0)) ** 0.5
#                         Amax_syn = np.max(np.sum(wpFu**2, axis=0)) ** 0.5
#                         Amax_syn1 = np.max(np.sum(wpFu1**2, axis=0)) ** 0.5
#                         win_b = win_starttime - syn_starttime
#                         win_e = win_endtime - syn_starttime
#                         obs_filt_proj = np.dot(proj_matrix, obs_filt)
#                         syn_filt_proj = np.dot(proj_matrix, syn_filt)
#                         for i in range(3):
#                             plt.subplot(311 + i)
#                             if i == 0:
#                                 title_str = "%s.%s " % (station_id, window_id)
#                                 title_str += "step_size:%.2f " % (step_size[istep])
#                                 title_str += "NCCwin:%.2f" % (cc_wpFd_wpFu1)
#                                 plt.title(title_str)
#
#                             idx_plt = range(syn_nl, (syn_nl + syn_npts))
#                             # whole trace, projected
#                             plt.plot(
#                                 syn_times[idx_plt],
#                                 obs_filt_proj[i, idx_plt] / Amax_obs,
#                                 "k",
#                                 linewidth=0.2,
#                             )
#                             plt.plot(
#                                 syn_times[idx_plt],
#                                 syn_filt_proj[i, idx_plt] / Amax_syn,
#                                 "b",
#                                 linewidth=0.2,
#                             )
#                             # windowed trace
#                             idx_plt = (win_b <= syn_times) & (syn_times <= win_e)
#                             plt.plot(
#                                 syn_times[idx_plt],
#                                 wpFd[i, idx_plt] / Amax_obs,
#                                 "k",
#                                 linewidth=1.0,
#                             )
#                             plt.plot(
#                                 syn_times[idx_plt],
#                                 wpFu[i, idx_plt] / Amax_syn,
#                                 "b",
#                                 linewidth=1.0,
#                             )
#                             plt.plot(
#                                 syn_times[idx_plt],
#                                 wpFu1[i, idx_plt] / Amax_syn1,
#                                 "r",
#                                 linewidth=1.0,
#                             )
#
#                             plt.xlim(
#                                 (syn_times[syn_nl], syn_times[syn_nl + syn_npts - 1])
#                             )
#                             plt.ylim((-1.5, 1.5))
#                             plt.ylabel(syn_orientation_codes[i])
#                         plt.show()
#             # end for window_id in window_dict:
#         # end for station_id in station_dict:
#
#         return wcc_sum, weight_sum
#

#     #
#     # ======================================================
#     #
#
#     def plot_seismogram_3comp(
#         self,
#         savefig=False,
#         out_dir="plot",
#         plot_param={
#             "time": [0, 100],
#             "rayp": 10.0,
#             "azbin": 10,
#             "window_id": "F.p,P",
#             "SNR": None,
#             "CC0": None,
#             "CCmax": None,
#             "dist": None,
#         },
#     ):
#         """Plot seismograms for one event
#         azbin: azimuthal bin size
#         win:
#         """
#         comp_name = ["R", "T", "Z"]
#         # ------ selection parameters
#         plot_time = plot_param["time"]
#         plot_azbin = plot_param["azbin"]
#         plot_rayp = plot_param["rayp"]
#         plot_window_id = plot_param["window_id"]
#
#         plot_SNR = np.array(plot_param["SNR"])
#         plot_CC0 = np.array(plot_param["CC0"])
#         plot_CCmax = np.array(plot_param["CCmax"])
#         plot_dist = np.array(plot_param["dist"])
#
#         # ------ event info
#         event = self.data["event"]
#         t0 = event["t0"]
#         tau = event["tau"]
#         evla = event["latitude"]
#         evlo = event["longitude"]
#         evdp = event["depth"]
#         mt = event["mt_rtp"]
#         Mrr = mt[0][0]
#         Mtt = mt[1][1]
#         Mpp = mt[2][2]
#         Mrt = mt[0][1]
#         Mrp = mt[0][2]
#         Mtp = mt[1][2]
#         focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
#
#         # ------ station info
#         station_dict = self.data["station"]
#         stla_all = []
#         stlo_all = []
#         dist_all = []
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             meta = station["meta"]
#             window_dict = station["window"]
#             # select data
#             if station["stat"]["code"] < 1:
#                 continue
#             if plot_window_id not in window_dict:
#                 continue
#             stla_all.append(meta["latitude"])
#             stlo_all.append(meta["longitude"])
#             dist_all.append(meta["dist_degree"])
#
#         # ------ traveltime curve
#         model = TauPyModel(model="ak135")
#         dist_ttcurve = np.arange(0.0, max(dist_all), 0.5)
#         ttcurve_p = []
#         ttcurve_P = []
#         ttcurve_s = []
#         ttcurve_S = []
#         for dist in dist_ttcurve:
#             arrivals = model.get_travel_times(
#                 source_depth_in_km=evdp,
#                 distance_in_degree=dist,
#                 phase_list=["p", "P", "s", "S"],
#             )
#             for arr in arrivals:
#                 if arr.name == "p":
#                     ttcurve_p.append((arr.distance, arr.time, arr.ray_param))
#                 elif arr.name == "P":
#                     ttcurve_P.append((arr.distance, arr.time, arr.ray_param))
#                 elif arr.name == "s":
#                     ttcurve_s.append((arr.distance, arr.time, arr.ray_param))
#                 elif arr.name == "S":
#                     ttcurve_S.append((arr.distance, arr.time, arr.ray_param))
#         # sort phases
#         ttcurve_p = sorted(ttcurve_p, key=lambda x: x[2])
#         ttcurve_P = sorted(ttcurve_P, key=lambda x: x[2])
#         ttcurve_s = sorted(ttcurve_s, key=lambda x: x[2])
#         ttcurve_S = sorted(ttcurve_S, key=lambda x: x[2])
#
#         # ------ map configuration
#         min_lat = min(min(stla_all), evla)
#         max_lat = max(max(stla_all), evla)
#         lat_range = max_lat - min_lat
#         min_lat -= 0.1 * lat_range
#         max_lat += 0.1 * lat_range
#         min_lon = min(min(stlo_all), evlo)
#         max_lon = max(max(stlo_all), evlo)
#         lon_range = max_lon - min_lon
#         min_lon -= 0.1 * lon_range
#         max_lon += 0.1 * lon_range
#         lat_0 = np.mean(stla_all)
#         lon_0 = np.mean(stlo_all)
#         #
#         parallels = np.arange(0.0, 81, 10.0)
#         meridians = np.arange(0.0, 351, 10.0)
#
#         # ------ plot azimuthal bins (one figure per azbin)
#         if plot_azbin <= 0.0:
#             raise Exception("plot_param['azbin']=%f must > 0.0" % plot_azbin)
#
#         for az in np.arange(0, 360, plot_azbin):
#             azmin = az
#             azmax = az + plot_azbin
#
#             print(azmin, azmax)
#
#             # ---- gather data for the current azbin
#             data_azbin = {}
#             for station_id in station_dict:
#                 # skip bad station
#                 station = station_dict[station_id]
#                 if station["stat"]["code"] < 1:
#                     continue
#
#                 # skip un-selected station
#                 meta = station["meta"]
#                 azimuth = meta["azimuth"]
#                 dist_degree = meta["dist_degree"]
#                 if plot_dist.any():
#                     if dist_degree < np.min(plot_dist) or dist_degree > np.max(
#                         plot_dist
#                     ):
#                         continue
#                 if azimuth < azmin or azimuth >= azmax:
#                     continue
#                 if plot_window_id not in window_dict:
#                     continue
#
#                 # skip bad window
#                 window_dict = station["window"]
#                 window = window_dict[plot_window_id]
#                 if window["stat"]["code"] <= 0:
#                     continue
#
#                 quality = window["quality"]
#                 if plot_SNR and quality["SNR"] < np.min(plot_SNR):
#                     continue
#
#                 cc = window["cc"]
#                 if plot_CC0 and cc["CC0"] < np.min(plot_CC0):
#                     continue
#                 if plot_CCmax and cc["CCmax"] < np.min(plot_CCmax):
#                     continue
#
#                 # get seismograms: syn/obs
#                 waveform = station["waveform"]
#                 time_sample = waveform["time_sample"]
#                 syn_starttime = time_sample["starttime"]
#                 syn_npts = time_sample["nt"]
#                 syn_delta = time_sample["delta"]
#                 syn_nyq = 0.5 / syn_delta
#                 obs = waveform["obs"]
#                 grn = waveform["grn"]
#                 # filter parameter
#                 filter_param = window["filter"]
#                 filter_a = filter_param["a"]
#                 filter_b = filter_param["b"]
#                 # filter seismograms
#                 obs = signal.filtfilt(filter_b, filter_a, obs)
#                 grn = signal.filtfilt(filter_b, filter_a, grn)
#                 # convolve stf on grn
#                 syn_freq = np.fft.rfftfreq(syn_npts, d=syn_delta)
#                 F_src = stf_gauss_spectrum(syn_freq, event["tau"])
#                 syn = np.fft.irfft(F_src * np.fft.rfft(grn), syn_npts)
#                 # project to polarity defined by the window
#                 proj_matrix = window["polarity"]["proj_matrix"]
#                 obs = np.dot(proj_matrix, obs)
#                 syn = np.dot(proj_matrix, syn)
#                 # rotate EN(0:2) -> RT(0:2) (T-R-Z: right-hand convention)
#                 Raz = (meta["back_azimuth"] + 180.0) % 360.0
#                 sin_Raz = np.sin(np.deg2rad(Raz))
#                 cos_Raz = np.cos(np.deg2rad(Raz))
#                 proj_matrix = [[sin_Raz, cos_Raz], [cos_Raz, -sin_Raz]]
#                 obs[0:2, :] = np.dot(proj_matrix, obs[0:2, :])
#                 syn[0:2, :] = np.dot(proj_matrix, syn[0:2, :])
#                 # append to data
#                 data_dict = {"meta": meta, "window": window, "syn": syn, "obs": obs}
#                 data_azbin[station_id] = data_dict
#             # endfor station_id in station_dict:
#
#             # ---- skip empty azbin
#             if not data_azbin:
#                 warn_str = "No station in the azbin [%f %f]." % (azmin, azmax)
#                 warnings.warn(warn_str)
#                 continue
#
#             # ---- create figure
#             fig = plt.figure(figsize=(8.5, 11))  # US letter
#             str_title = "{:s} (win: {:s}, az: {:04.1f}~{:04.1f})".format(
#                 event["id"], plot_window_id, azmin, azmax
#             )
#             fig.text(0.5, 0.95, str_title, size="x-large", horizontalalignment="center")
#             # ---- plot station/event map
#             ax_origin = [0.3, 0.74]
#             ax_size = [0.4, 0.2]
#             ax_map = fig.add_axes(ax_origin + ax_size)
#             ax_bm = Basemap(
#                 projection="merc",
#                 resolution="l",
#                 llcrnrlat=min_lat,
#                 llcrnrlon=min_lon,
#                 urcrnrlat=max_lat,
#                 urcrnrlon=max_lon,
#                 lat_0=lat_0,
#                 lon_0=lon_0,
#             )
#             ax_bm.drawcoastlines(linewidth=0.1)
#             ax_bm.drawcountries(linewidth=0.1)
#             ax_bm.drawparallels(
#                 parallels, linewidth=0.1, labels=[1, 0, 0, 1], fontsize=10, fmt="%3.0f"
#             )
#             ax_bm.drawmeridians(
#                 meridians, linewidth=0.1, labels=[1, 0, 0, 1], fontsize=10, fmt="%3.0f"
#             )
#             sx, sy = ax_bm(stlo_all, stla_all)
#             ax_bm.scatter(sx, sy, s=10, marker="^", facecolor="blue", edgecolor="")
#             # plot focal mechanism
#             sx, sy = ax_bm(evlo, evla)
#             bb_width = 110000.0 * np.abs(max(stlo_all) - min(stlo_all)) * 0.1
#             b = Beach(focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor="r")
#             ax_map.add_collection(b)
#             # -- plot the station location
#             stla = [x["meta"]["latitude"] for x in data_azbin.itervalues()]
#             stlo = [x["meta"]["longitude"] for x in data_azbin.itervalues()]
#             sx, sy = ax_bm(stlo, stla)
#             ax_bm.scatter(sx, sy, s=10, marker="^", facecolor="red", edgecolor="")
#
#             # -- create axis for seismograms
#             ax_RTZ = []
#             for i in range(3):
#                 ax_origin = [0.07 + 0.3 * i, 0.05]
#                 ax_size = [0.25, 0.65]
#                 ax_RTZ.append(fig.add_axes(ax_origin + ax_size))
#
#             # -- plot traveltime curves
#             for i in range(3):
#                 ax = ax_RTZ[i]
#                 ax.plot(
#                     [x[1] - plot_rayp * x[0] for x in ttcurve_p],
#                     [x[0] for x in ttcurve_p],
#                     "b-",
#                     linewidth=0.2,
#                 )
#                 ax.plot(
#                     [x[1] - plot_rayp * x[0] for x in ttcurve_P],
#                     [x[0] for x in ttcurve_P],
#                     "b-",
#                     linewidth=0.2,
#                 )
#                 ax.plot(
#                     [x[1] - plot_rayp * x[0] for x in ttcurve_s],
#                     [x[0] for x in ttcurve_s],
#                     "c-",
#                     linewidth=0.2,
#                 )
#                 ax.plot(
#                     [x[1] - plot_rayp * x[0] for x in ttcurve_S],
#                     [x[0] for x in ttcurve_S],
#                     "c-",
#                     linewidth=0.2,
#                 )
#
#             # -- ylim setting
#             y = [x["meta"]["dist_degree"] for x in data_azbin.itervalues()]
#             ny = len(y)
#             plot_dy = 0.5 * (max(y) - min(y) + 1) / ny
#             if plot_dist:
#                 plot_ymax = max(plot_dist) + 2 * plot_dy
#                 plot_ymin = min(plot_dist) - 2 * plot_dy
#             else:
#                 plot_ymax = max(y) + 2 * plot_dy
#                 plot_ymin = min(y) - 2 * plot_dy
#
#             # -- plot each station
#             for station_id in data_azbin:
#                 station = data_azbin[station_id]
#                 meta = station["meta"]
#                 window = station["window"]
#                 syn = station["syn"]
#                 obs = station["obs"]
#
#                 # get plot time
#                 dist_degree = meta["dist_degree"]
#                 reduced_time = dist_degree * plot_rayp
#                 # time of first sample referred to centroid time
#                 t0 = syn_starttime - event["t0"]
#                 # time of samples referred to centroid time
#                 syn_times = syn_delta * np.arange(syn_npts) + t0
#                 # plot time window
#                 plot_t0 = min(plot_time) + reduced_time
#                 plot_t1 = max(plot_time) + reduced_time
#                 plot_idx = (syn_times > plot_t0) & (syn_times < plot_t1)
#                 # plot time
#                 t_plot = syn_times[plot_idx] - reduced_time
#
#                 #  window begin/end
#                 taper = window["taper"]
#                 win_starttime = taper["starttime"]
#                 win_endtime = taper["endtime"]
#                 win_t0 = (win_starttime - event["t0"]) - reduced_time
#                 win_t1 = (win_endtime - event["t0"]) - reduced_time
#
#                 # plot seismograms
#                 Amax_obs = np.sqrt(np.max(np.sum(obs[:, plot_idx] ** 2, axis=0)))
#                 Amax_syn = np.sqrt(np.max(np.sum(syn[:, plot_idx] ** 2, axis=0)))
#                 for i in range(3):
#                     ax = ax_RTZ[i]
#                     ax.plot(
#                         t_plot,
#                         plot_dy * obs[i, plot_idx] / Amax_obs + dist_degree,
#                         "k-",
#                         linewidth=0.5,
#                     )
#                     ax.plot(
#                         t_plot,
#                         plot_dy * syn[i, plot_idx] / Amax_syn + dist_degree,
#                         "r-",
#                         linewidth=0.5,
#                     )
#                     # mark measure window range
#                     ax.plot(win_t0, dist_degree, "k|", markersize=8)
#                     ax.plot(win_t1, dist_degree, "k|", markersize=8)
#                     # annotate amplitude
#                     if i == 0:
#                         ax.text(
#                             max(plot_time),
#                             dist_degree,
#                             "%.1e " % (Amax_obs),
#                             verticalalignment="bottom",
#                             horizontalalignment="right",
#                             fontsize=7,
#                             color="black",
#                         )
#                         ax.text(
#                             max(plot_time),
#                             dist_degree,
#                             "%.1e " % (Amax_syn),
#                             verticalalignment="top",
#                             horizontalalignment="right",
#                             fontsize=7,
#                             color="red",
#                         )
#                     # annotate CC0
#                     if i == 0:
#                         ax.text(
#                             max(plot_time),
#                             dist_degree,
#                             " %.3f" % (window["cc"]["CC0"]),
#                             verticalalignment="center",
#                             fontsize=7,
#                         )
#                     # annotate window weight
#                     if i == 1:
#                         ax.text(
#                             max(plot_time),
#                             dist_degree,
#                             " %.1f" % (window["weight"]),
#                             verticalalignment="center",
#                             fontsize=7,
#                         )
#                     # annotate station names
#                     if i == 2:
#                         # str_annot = '%.3f,%.1f,%s' % (
#                         #    misfit['CC0'], window['weight'], station_id)
#                         ax.text(
#                             max(plot_time),
#                             dist_degree,
#                             " " + station_id,
#                             verticalalignment="center",
#                             fontsize=7,
#                         )
#             # endfor data in data_azbin:
#
#             # -- set axes limits and lables, annotation
#             for i in range(3):
#                 ax = ax_RTZ[i]
#                 ax.set_xlim(min(plot_time), max(plot_time))
#                 ax.set_ylim(plot_ymin, plot_ymax)
#                 ax.set_title(comp_name[i])
#                 ax.set_xlabel("t - {:.1f}*dist (s)".format(plot_rayp))
#                 ax.tick_params(axis="both", labelsize=10)
#                 # ylabel
#                 if i == 0:
#                     ax.set_ylabel("dist (deg)")
#                 else:
#                     ax.set_yticklabels([])
#
#             # -- save figures
#             if savefig:
#                 out_file = "%s/%s_az_%03d_%03d_%s.pdf" % (
#                     out_dir,
#                     event["id"],
#                     azmin,
#                     azmax,
#                     plot_window_id,
#                 )
#                 plt.savefig(out_file, format="pdf")
#             else:
#                 plt.show()
#             plt.close(fig)
#
#     #
#     # ======================================================
#     #
#
#     def plot_seismogram_1comp(
#         self,
#         savefig=False,
#         out_dir="plot",
#         window_id="p,P_Z",
#         azbin=10,
#         begin_time=0,
#         end_time=0,
#         clip_ratio=1.5,
#         min_CC0=None,
#         min_CCmax=None,
#         min_SNR=None,
#         dist_lim=None,
#         plot_az0=0,
#         plot_adj=False,  # whether plot adjoint source
#         align_time=False,  # whether align the phase according to cc time shift
#     ):
#         """
#         Plot record section in azimuthal bins
#
#         Parameters
#         ----------
#         azbin: azimuthal bin size
#
#         begin/end_time: time range that is added to the automatically determined
#           plot time range. See below.
#
#         clip: do not plot waveform with amplitudes larger than
#           <clip>*max_amplitude_in_select_time_window
#
#         Notes
#         -----
#           The time for plot is reduced time relative to origin time + rayp*dist
#
#           linear regression is done to find the average rayp of misfit windows
#           and the plot begin/end time is found that can include all misfit windows
#           in the plot.
#
#         """
#         # ------ check parameters
#         # plot_time = np.array([begin_time, end_time])
#
#         # in case reverse the distance axis
#         # plot_flip = -1
#         plot_flip = 1
#
#         plot_azbin = float(azbin)
#         if plot_azbin <= 0:
#             raise Exception("plot_azbin(%f) should be larger than 0.0" % (plot_azbin))
#
#         plot_window_id = window_id
#         plot_SNR = np.array(min_SNR)
#         plot_CC0 = np.array(min_CC0)
#         plot_CCmax = np.array(min_CCmax)
#         plot_dist = np.array(dist_lim)
#
#         plot_clip = float(clip_ratio)
#         if plot_clip < 1.0:
#             raise Exception("clip_ratio(%f) should be larger than 1.0" % (plot_clip))
#
#         # ------ event info
#         event = self.data["event"]
#         t0 = event["t0"]
#         tau = event["tau"]
#         evla = event["latitude"]
#         evlo = event["longitude"]
#         evdp = event["depth"]
#         # evdp has to be >=0 otherwise taup would crash
#         if evdp < 0.0:
#             evdp = 0.0
#
#         mt = event["mt_rtp"]
#         Mrr = mt[0][0]
#         Mtt = mt[1][1]
#         Mpp = mt[2][2]
#         Mrt = mt[0][1]
#         Mrp = mt[0][2]
#         Mtp = mt[1][2]
#         focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
#
#         # ------ get station info
#         station_dict = self.data["station"]
#         stla_all = []
#         stlo_all = []
#         dist_all = []
#         winb_all = []
#         wine_all = []
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             meta = station["meta"]
#             window_dict = station["window"]
#             # select data
#             if station["stat"]["code"] < 1:
#                 continue
#             if plot_window_id not in window_dict:
#                 continue
#             stla_all.append(meta["latitude"])
#             stlo_all.append(meta["longitude"])
#             dist_all.append(meta["dist_degree"])
#             taper = window_dict[plot_window_id]["taper"]
#             winb_all.append(taper["starttime"] - event["t0"])
#             wine_all.append(taper["endtime"] - event["t0"])
#
#         if not dist_all:
#             warnings.warn("No data to plot!")
#             return
#
#         # get average moveout of the window center
#         dist_all = np.array(dist_all)
#         winb_all = np.array(winb_all)
#         wine_all = np.array(wine_all)
#         winc_all = (winb_all + wine_all) / 2.0
#         # linear regression tc = dist*rayp + tb
#         A = np.vstack([dist_all, np.ones(len(dist_all))]).T
#         plot_rayp, plot_c = np.linalg.lstsq(A, winc_all)[0]
#         # round to the integer
#         plot_rayp = np.round(plot_rayp)
#         # plot_rayp = 16 # taok: temporary use
#         # KT KT this should be decided for each azimuthal bin. So I moved
#         #       the following 6 lines into the "plot azimuthal bin" section.
#         ## get time window relative to the regressed window central time
#         # plot_t0 = np.min(winb_all - plot_rayp*dist_all)
#         # plot_t1 = np.max(wine_all - plot_rayp*dist_all)
#         ## modify the plot time rage
#         # plot_time[0] += plot_t0
#         # plot_time[1] += plot_t1
#
#         # ------ calculate traveltime curves (only for body wave)
#         phase_names = plot_window_id.split("_")[0]
#         if phase_names not in ["surface", "Rayleigh", "Love"]:
#             model = TauPyModel(model="ak135")
#             # distance samples
#             dist_ttcurve = np.arange(0.0, max(dist_all), 0.5)
#             phase_list = [x for x in phase_names.split(",")]
#             ttcurve = {}
#             for phase_name in phase_list:
#                 ttcurve[phase_name] = []
#             for dist in dist_ttcurve:
#                 arrivals = model.get_travel_times(
#                     source_depth_in_km=evdp,
#                     distance_in_degree=dist,
#                     phase_list=phase_list,
#                 )
#                 for arr in arrivals:
#                     for phase_name in phase_list:
#                         if arr.name == phase_name:
#                             ttcurve[phase_name].append(
#                                 (arr.distance, arr.time, arr.ray_param)
#                             )
#             # sort (dist, ttime, rayp) points based on ray parameter
#             for phase_name in phase_list:
#                 ttcurve[phase_name] = sorted(ttcurve[phase_name], key=lambda x: x[2])
#
#         # ------ map configuration
#         min_lat = min(min(stla_all), evla)
#         max_lat = max(max(stla_all), evla)
#         lat_range = max_lat - min_lat
#         min_lat -= 0.1 * lat_range
#         max_lat += 0.1 * lat_range
#         if min_lat < -90.0:
#             min_lat = -90.0
#         if max_lat > 90.0:
#             max_lat = 90.0
#         min_lon = min(min(stlo_all), evlo)
#         max_lon = max(max(stlo_all), evlo)
#         lon_range = max_lon - min_lon
#         min_lon -= 0.1 * lon_range
#         max_lon += 0.1 * lon_range
#         if min_lon < -180.0:
#             min_lon = -180.0
#         if max_lon > 180.0:
#             max_lon = 180.0
#         lat_0 = np.mean(stla_all)
#         lon_0 = np.mean(stlo_all)
#         #
#         parallels = np.arange(-90.0, 89.0, 10.0)
#         meridians = np.arange(0.0, 351, 10.0)
#
#         # ------ plot azimuthal bins (one figure per azbin)
#         if plot_azbin <= 0.0:
#             raise Exception("plot_param['azbin']=%f must > 0.0" % plot_azbin)
#
#         for az in np.arange(plot_az0, plot_az0 + 360, plot_azbin):
#             azmin = az
#             azmax = az + plot_azbin
#
#             print("Azimuthal range: ", azmin, azmax)
#
#             # ---- gather data for the current azbin
#             data_azbin = {}
#             for station_id in station_dict:
#                 station = station_dict[station_id]
#                 # skip bad station
#                 if station["stat"]["code"] < 1:
#                     continue
#                 # skip station not in the selection criteria
#                 meta = station["meta"]
#                 azimuth = meta["azimuth"]
#                 dist_degree = meta["dist_degree"]
#                 if plot_dist.any():
#                     if dist_degree < np.min(plot_dist) or dist_degree > np.max(
#                         plot_dist
#                     ):
#                         continue
#                 # if azimuth < azmin or azimuth >= azmax:
#                 if (azimuth - azmin) % 360 >= plot_azbin:
#                     continue
#
#                 window_dict = station["window"]
#                 # check if required window exists
#                 if plot_window_id not in window_dict:
#                     continue
#                 window = window_dict[plot_window_id]
#                 # skip bad window
#                 if window["stat"]["code"] <= 0:
#                     continue
#                 # skip window which does not pass the selection criteria
#                 quality = window["quality"]
#                 if plot_SNR and quality["SNR"] < np.min(plot_SNR):
#                     continue
#                 cc = window["cc"]
#                 if plot_CC0 and cc["CC0"] < np.min(plot_CC0):
#                     continue
#                 if plot_CCmax and cc["CCmax"] < np.min(plot_CCmax):
#                     continue
#
#                 # get seismograms: syn/obs
#                 waveform = station["waveform"]
#                 time_sample = waveform["time_sample"]
#                 syn_starttime = time_sample["starttime"]
#                 syn_npts = time_sample["nt"]
#                 syn_delta = time_sample["delta"]
#                 syn_nyq = 0.5 / syn_delta
#                 # filter parameter
#                 filter_param = window["filter"]
#                 filter_a = filter_param["a"]
#                 filter_b = filter_param["b"]
#                 # filter seismograms
#                 obs = signal.filtfilt(filter_b, filter_a, waveform["obs"])
#                 if "syn" in waveform:
#                     syn = signal.filtfilt(filter_b, filter_a, waveform["syn"])
#                 elif "grn" in waveform:
#                     grn = signal.filtfilt(filter_b, filter_a, waveform["grn"])
#                     # convolve stf with grn
#                     syn_freq = np.fft.rfftfreq(syn_npts, d=syn_delta)
#                     F_src = stf_gauss_spectrum(syn_freq, event["tau"])
#                     syn = np.fft.irfft(F_src * np.fft.rfft(grn), syn_npts)
#                 else:
#                     err = "station(%s) has no syn or grn in waveform data." % (
#                         station_id
#                     )
#                     raise Exception(err)
#                 # project to polarity defined by the window
#                 polarity = window["polarity"]
#                 comp = polarity["component"]
#                 cmpaz = polarity["azimuth"]
#                 cmpdip = polarity["dip"]
#                 if comp in ["Z", "R", "T"]:
#                     sin_az = np.sin(np.deg2rad(cmpaz))
#                     cos_az = np.cos(np.deg2rad(cmpaz))
#                     sin_dip = np.sin(np.deg2rad(cmpdip))
#                     cos_dip = np.cos(np.deg2rad(cmpdip))
#                     cmp_vec = np.array(
#                         [
#                             cos_dip * sin_az,  # cos(E, comp)
#                             cos_dip * cos_az,  # N, comp
#                             -sin_dip,
#                         ]
#                     )  # Z, comp
#                 else:
#                     raise Exception("Not single component: " % (comp))
#                 obs = np.dot(cmp_vec, obs)
#                 syn = np.dot(cmp_vec, syn)
#
#                 # append to data
#                 if plot_adj:
#                     adj = np.dot(cmp_vec, station["dchi_du"])
#                     data_dict = {
#                         "meta": meta,
#                         "window": window,
#                         "syn": syn,
#                         "obs": obs,
#                         "adj": adj,
#                     }
#                 else:
#                     data_dict = {
#                         "meta": meta,
#                         "window": window,
#                         "syn": syn,
#                         "obs": obs,
#                     }
#                 data_azbin[station_id] = data_dict
#             # endfor station_id in station_dict:
#
#             # ---- skip empty azbin
#             if not data_azbin:
#                 warn_str = "No station in the azbin [%f %f]." % (azmin, azmax)
#                 warnings.warn(warn_str)
#                 continue
#
#             # ---- create figure
#             fig = plt.figure(figsize=(11, 8.5))  # US Letter
#             str_title = "{:s} ({:s} az:{:04.1f}~{:04.1f} dep:{:.1f})".format(
#                 event["id"], plot_window_id, azmin, azmax, event["depth"]
#             )
#             fig.text(
#                 0.5, 0.965, str_title, size="x-large", horizontalalignment="center"
#             )
#
#             # ---- plot station/event map
#             ax_origin = [0.05, 0.60]
#             ax_size = [0.3, 0.3]
#             ax_map = fig.add_axes(ax_origin + ax_size)
#             ax_bm = Basemap(
#                 projection="merc",
#                 resolution="l",
#                 llcrnrlat=min_lat,
#                 llcrnrlon=min_lon,
#                 urcrnrlat=max_lat,
#                 urcrnrlon=max_lon,
#                 lat_0=lat_0,
#                 lon_0=lon_0,
#             )
#             ax_bm.drawcoastlines(linewidth=0.1)
#             ax_bm.drawcountries(linewidth=0.1)
#             ax_bm.drawparallels(
#                 parallels, linewidth=0.1, labels=[1, 0, 0, 1], fontsize=10, fmt="%3.0f"
#             )
#             ax_bm.drawmeridians(
#                 meridians, linewidth=0.1, labels=[1, 0, 0, 1], fontsize=10, fmt="%3.0f"
#             )
#             sx, sy = ax_bm(stlo_all, stla_all)
#             ax_bm.scatter(sx, sy, s=10, marker="^", facecolor="blue", edgecolor="")
#             # plot focal mechanism
#             sx, sy = ax_bm(evlo, evla)
#             bb_width = 110000.0 * np.abs(max(stlo_all) - min(stlo_all)) * 0.1
#             b = beach(focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor="r")
#             ax_map.add_collection(b)
#             # -- plot the station location
#             stla = [data_azbin[key]["meta"]["latitude"] for key in data_azbin]
#             stlo = [data_azbin[key]["meta"]["longitude"] for key in data_azbin]
#             sx, sy = ax_bm(stlo, stla)
#             ax_bm.scatter(sx, sy, s=10, marker="^", facecolor="red", edgecolor="")
#
#             # -- create axis for seismograms
#             ax_origin = [0.45, 0.05]
#             ax_size = [0.43, 0.90]
#             # ax_size = [0.3, 0.90]
#             ax_1comp = fig.add_axes(ax_origin + ax_size)
#
#             # -- xlim setting
#             win_all = [data_azbin[key]["window"] for key in data_azbin]
#             winb_all = np.array(
#                 [win["taper"]["starttime"] - event["t0"] for win in win_all]
#             )
#             wine_all = np.array(
#                 [win["taper"]["endtime"] - event["t0"] for win in win_all]
#             )
#             dist_all = np.array(
#                 [data_azbin[key]["meta"]["dist_degree"] for key in data_azbin]
#             )
#             # get time window relative to the regressed window central time
#             plot_t0 = np.min(winb_all - plot_rayp * dist_all)
#             plot_t1 = np.max(wine_all - plot_rayp * dist_all)
#             plot_time = np.array([begin_time + plot_t0, end_time + plot_t1])
#
#             # -- ylim setting
#             y = [data_azbin[key]["meta"]["dist_degree"] for key in data_azbin]
#             ny = len(y)
#             plot_dy = 0.5 * (max(y) - min(y) + 1) / ny
#             if plot_dist.any():
#                 plot_ymax = max(plot_dist) + 2 * plot_dy
#                 plot_ymin = min(plot_dist) - 2 * plot_dy
#             else:
#                 plot_ymax = max(y) + 2 * plot_dy
#                 plot_ymin = min(y) - 2 * plot_dy
#
#             # -- plot traveltime curves
#             if phase_names not in ["surface", "Rayleigh", "Love"]:
#                 for phase_name in phase_list:
#                     # skip if no tt curves for this phase_names
#                     if not ttcurve[phase_name]:
#                         continue
#                     # reduced time
#                     phase_times = np.array(
#                         [x[1] - plot_rayp * x[0] for x in ttcurve[phase_name]]
#                     )
#                     phase_distances = np.array([x[0] for x in ttcurve[phase_name]])
#                     # skip if not in plot range
#                     max_dist = np.max(phase_distances)
#                     min_dist = np.min(phase_distances)
#                     if max_dist < plot_ymin or min_dist > plot_ymax:
#                         continue
#                     ax_1comp.plot(phase_times, phase_distances, "b-", linewidth=0.1)
#                     # ax_1comp.plot(phase_times, phase_distances, 'b.', markersize=0.5)
#                     # label phase names
#                     if max_dist < plot_ymax:
#                         y_str = max_dist
#                         x_str = max(phase_times[phase_distances == max_dist])
#                     else:
#                         y_str = plot_ymax
#                         max_dist = max(phase_distances[phase_distances <= plot_ymax])
#                         x_str = max(phase_times[phase_distances == max_dist])
#                     ax_1comp.text(
#                         x_str,
#                         y_str,
#                         phase_name,
#                         verticalalignment="top",
#                         horizontalalignment="center",
#                         fontsize=11,
#                         color="blue",
#                     )
#
#             # -- plot each station
#             if plot_adj:  # use a constant scaling factor for adj_src
#                 Amax_adj = -1.0
#                 for station_id in data_azbin:
#                     station = data_azbin[station_id]
#                     meta = station["meta"]
#                     window = station["window"]
#                     adj = station["adj"]
#                     # get plot time
#                     dist_degree = meta["dist_degree"]
#                     reduced_time = dist_degree * plot_rayp
#                     # time of first sample referred to centroid time
#                     t0 = syn_starttime - event["t0"]
#                     # time of samples referred to centroid time
#                     syn_times = syn_delta * np.arange(syn_npts) + t0
#                     ## plot time window
#                     # plot_t0 = min(plot_time) + reduced_time
#                     # plot_t1 = max(plot_time) + reduced_time
#                     # plot_idx = (syn_times > plot_t0) & (syn_times < plot_t1)
#                     ## plot time (reduced time)
#                     # t_plot = syn_times[plot_idx] - reduced_time
#                     #  window begin/end
#                     taper = window["taper"]
#                     win_starttime = taper["starttime"] - event["t0"]
#                     win_endtime = taper["endtime"] - event["t0"]
#                     win_t0 = win_starttime - reduced_time
#                     win_t1 = win_endtime - reduced_time
#                     win_idx = (syn_times > win_starttime) & (syn_times < win_endtime)
#
#                     Amax_adj = max(Amax_adj, np.sqrt(np.max(adj[win_idx] ** 2)))
#
#             for station_id in data_azbin:
#                 station = data_azbin[station_id]
#                 meta = station["meta"]
#                 window = station["window"]
#                 syn = station["syn"]
#                 obs = station["obs"]
#
#                 if align_time:
#                     cc_tshift = window["cc"]["cc_tshift"]
#
#                 # get plot time
#                 dist_degree = meta["dist_degree"]
#                 reduced_time = dist_degree * plot_rayp
#                 # time of first sample referred to centroid time
#                 t0 = syn_starttime - event["t0"]
#                 # time of samples referred to centroid time
#                 syn_times = syn_delta * np.arange(syn_npts) + t0
#                 # plot time window
#                 plot_t0 = min(plot_time) + reduced_time
#                 plot_t1 = max(plot_time) + reduced_time
#                 plot_idx = (syn_times > plot_t0) & (syn_times < plot_t1)
#                 # plot time (reduced time)
#                 t_plot = syn_times[plot_idx] - reduced_time
#
#                 #  window begin/end
#                 taper = window["taper"]
#                 win_starttime = taper["starttime"] - event["t0"]
#                 win_endtime = taper["endtime"] - event["t0"]
#                 win_t0 = win_starttime - reduced_time
#                 win_t1 = win_endtime - reduced_time
#                 win_idx = (syn_times > win_starttime) & (syn_times < win_endtime)
#
#                 # plot seismograms
#                 Amax_obs = np.sqrt(np.max(obs[win_idx] ** 2))
#                 Amax_syn = np.sqrt(np.max(syn[win_idx] ** 2))
#
#                 # clip large amplitudes
#                 if plot_adj:
#                     adj = station["adj"]
#                     y = adj[plot_idx] / Amax_adj
#                     idx = abs(y) > plot_clip + 1.0e-3
#                     y[idx] = np.nan
#                     ax_1comp.plot(
#                         t_plot,
#                         plot_flip * plot_dy * y + dist_degree,
#                         "k-",
#                         linewidth=0.5,
#                     )
#
#                 y = obs[plot_idx] / Amax_obs
#                 idx = abs(y) > plot_clip + 1.0e-3
#                 y[idx] = np.nan
#                 ax_1comp.plot(
#                     t_plot, plot_flip * plot_dy * y + dist_degree, "k-", linewidth=0.5
#                 )
#
#                 y = syn[plot_idx] / Amax_syn
#                 idx = abs(y) > plot_clip + 1.0e-3
#                 y[idx] = np.nan
#                 if align_time:
#                     ax_1comp.plot(
#                         t_plot + cc_tshift,
#                         plot_flip * plot_dy * y + dist_degree,
#                         "r-",
#                         linewidth=0.5,
#                     )
#                 else:
#                     ax_1comp.plot(
#                         t_plot,
#                         plot_flip * plot_dy * y + dist_degree,
#                         "r-",
#                         linewidth=0.5,
#                     )
#
#                 # mark measure window range
#                 ax_1comp.plot(win_t0, dist_degree, "k|", markersize=8)
#                 ax_1comp.plot(win_t1, dist_degree, "k|", markersize=8)
#                 ## annotate amplitude
#                 #  ax.text(max(plot_time), dist_degree, '%.1e ' % (Amax_obs),
#                 #      verticalalignment='bottom',
#                 #      horizontalalignment='right',
#                 #      fontsize=7, color='black')
#                 #  ax.text(max(plot_time), dist_degree, '%.1e ' % (Amax_syn),
#                 #      verticalalignment='top',
#                 #      horizontalalignment='right',
#                 #      fontsize=7, color='red')
#                 ## annotate CC0
#                 #  ax.text(max(plot_time), dist_degree, ' %.3f'%(window['cc']['CC0']),
#                 #      verticalalignment='center', fontsize=7)
#                 ## annotate window weight
#                 # if i == 1:
#                 #  ax.text(max(plot_time), dist_degree, ' %.1f' % (window['weight']),
#                 #      verticalalignment='center', fontsize=7)
#                 ##annotate station names
#                 str_annot = " %s (%.3f,%.3f,%.1f)" % (
#                     station_id,
#                     window["cc"]["CC0"],
#                     window["cc"]["cc_tshift"],
#                     window["weight"],
#                 )
#                 ax_1comp.text(
#                     max(plot_time),
#                     dist_degree,
#                     str_annot,
#                     verticalalignment="center",
#                     fontsize=7,
#                 )
#                 # ax_1comp.text(160, dist_degree, str_annot,
#                 #    verticalalignment='center', fontsize=7)
#
#             # endfor data in data_azbin:
#
#             # -- set axes limits and lables, annotation
#             ax_1comp.set_xlim(min(plot_time), max(plot_time))
#             # ax_1comp.set_xlim(80,160)
#             ax_1comp.set_ylim(plot_ymin, plot_ymax)
#             ax_1comp.set_xlabel("t - {:.1f}*dist (s)".format(plot_rayp))
#             ax_1comp.tick_params(axis="both", labelsize=10)
#             # ylabel
#             ax_1comp.set_ylabel("dist (deg)")
#             # ax_1comp.invert_yaxis()
#
#             # -- save figures
#             if savefig:
#                 out_file = "%s/%s_az_%03d_%03d_%s.pdf" % (
#                     out_dir,
#                     event["id"],
#                     azmin,
#                     azmax,
#                     plot_window_id,
#                 )
#                 plt.savefig(out_file, format="pdf")
#             else:
#                 plt.show()
#             plt.close(fig)
#
#     #
#     # ======================================================
#     #
#
#     def read_perturbed_waveform(
#         self,
#         syn_dir="output_perturb/sac",
#         syn_band_code="MX",
#         syn_suffix=".sem.sac",
#         model_name="perturb",
#         sac_dir=None,
#     ):
#         """
#         read in perturbed seismograms
#
#         """
#         syn_orientation_codes = ["E", "N", "Z"]
#
#         event = self.data["event"]
#         t0 = event["t0"]
#
#         station_dict = self.data["station"]
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected statations
#             if station["stat"]["code"] < 1:
#                 continue
#
#             # ------ time samples
#             waveform = station["waveform"]
#
#             if "syn" not in waveform:
#                 raise Exception(
#                     "%s: initial syn is not stored in waveform!" % (station_id)
#                 )
#
#             time_sample = waveform["time_sample"]
#             starttime = time_sample["starttime"]
#             dt = time_sample["delta"]
#             nt = time_sample["nt"]
#             nl = time_sample["nl"]  # npts of left padding
#             nr = time_sample["nr"]  # npts of right padding
#             sem_nt = nt - nl - nr  # number of time sample number in SEM simulation
#             t = np.arange(nt) * dt + (starttime - t0)  # referred to t0
#
#             # ------ get file paths of syn seismograms
#             syn_files = [
#                 "{:s}/{:s}.{:2s}{:1s}{:s}".format(
#                     syn_dir, station_id, syn_band_code, x, syn_suffix
#                 )
#                 for x in syn_orientation_codes
#             ]
#
#             # ------ read in syn seismograms
#             syn_st = read(syn_files[0])
#             syn_st += read(syn_files[1])
#             syn_st += read(syn_files[2])
#
#             # ------ check the same time samples as original syn
#             if not is_equal(
#                 [(tr.stats.starttime, tr.stats.delta, tr.stats.npts) for tr in syn_st]
#             ):
#                 raise Exception(
#                     "%s: not equal time samples in"
#                     " synthetic seismograms." % (station_id)
#                 )
#             tr = syn_st[0]
#
#             if tr.stats.delta != dt:
#                 raise Exception("%s: not the same dt for diff-srcloc!" % (station_id))
#
#             tr_starttime = tr.stats.starttime - nl * dt
#             if tr_starttime != starttime:
#                 raise Exception(
#                     "%s: not the same starttime for diff-srcloc!" % (station_id)
#                 )
#
#             if tr.stats.npts != sem_nt:
#                 raise Exception("%s: not the same npts for diff-srcloc!" % (station_id))
#
#             # ------ store perturbed waveform
#             syn_ENZ = np.zeros((3, nt))
#             for i in range(3):
#                 syn_ENZ[i, nl : (nl + sem_nt)] = syn_st[i].data
#             waveform[model_name] = syn_ENZ
#
#             # DEBUG: check du
#             # print(dchi
#             # for i in range(3):
#             #  plt.subplot(311+i)
#             #  #plt.plot(t,grn0[i,:],'k', t,syn_ENZ[i,:],'r', t,dg[i,:], 'b')
#             #  plt.plot(t, du[i,:], 'k')
#             # plt.show()
#             # if sac_dir:
#             #  for i in range(3):
#             #    tr.data = syn_ENZ[i,:]
#             #    tr.stats.starttime = starttime
#             #    tr.stats.delta = dt
#             #    tr.stats.npts = nt
#             #    out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
#             #        sac_dir, station_id, syn_band_code,
#             #        syn_orientation_codes[i])
#             #    tr.write(out_file, 'sac')
#
#     #
#     # ======================================================
#     #
#
#     def output_adj_for_perturbed_waveform(
#         self,
#         plot=False,
#         out_dir="adj",
#         model_name="perturb",
#         syn_band_code="MX",
#     ):
#         """
#         calculate adjoint sources (dchi_du) for perturbed waveform
#
#         Notes
#         -----
#         chi : misfit functional (normalized zero-lag correlation coef.)
#         u : synthetic waveform
#
#         """
#         # ------
#         syn_orientation_codes = ["E", "N", "Z"]
#
#         event = self.data["event"]
#         station_dict = self.data["station"]
#
#         tr = Trace()
#
#         # ====== loop each station
#         for station_id in station_dict:
#             station = station_dict[station_id]
#             # skip rejected statations
#             if station["stat"]["code"] < 1:
#                 continue
#
#             meta = station["meta"]
#             window_dict = station["window"]
#
#             waveform = station["waveform"]
#             time_sample = waveform["time_sample"]
#             syn_starttime = time_sample["starttime"]
#             syn_delta = time_sample["delta"]
#             syn_nt = time_sample["nt"]
#             syn_nl = time_sample["nl"]
#             syn_nr = time_sample["nr"]
#
#             # time samples for ascii output, referred to origin time
#             # without padding
#             npts = syn_nt - syn_nl - syn_nr
#             starttime = syn_starttime + syn_nl * syn_delta
#             syn_times = np.arange(npts) * syn_delta
#             b = starttime - event["t0"]
#             syn_times += b
#
#             # get obs and perturbed syn
#             obs = waveform["obs"]
#             if model_name in waveform:
#                 syn = waveform[model_name]
#             else:
#                 error_str = "%s: No synthetics found for %s" % (station_id, model_name)
#                 raise Exception(error_str)
#
#             # ------ loop each window
#             # sum of adjoint sources from all windows
#             dchi_du = np.zeros((3, syn_nt))
#
#             for window_id in window_dict:
#                 # window parameters
#                 window = window_dict[window_id]
#                 # skip windows not used in adjoint source from the unperturbed waveform
#                 if window["stat"]["code"] != 1:
#                     continue
#
#                 # ------ window parameters
#                 # filter
#                 filter_dict = window["filter"]
#                 filter_a = filter_dict["a"]
#                 filter_b = filter_dict["b"]
#                 # taper
#                 win_func = window["taper"]["win"]
#                 # polarity projection
#                 proj_matrix = window["polarity"]["proj_matrix"]
#                 # window weight
#                 weight = window["weight"]
#                 # misfit type
#                 misfit_type = window["misfit_type"]
#
#                 # ------ filter obs, syn
#                 # NOTE: use lfilter (causal filter) to avoid contamination from the right
#                 # end of the signal, but with asymmetric response and
#                 # peak shift ~ 1/4 min. period (e.g. 0.01-0.1Hz -> 2.5s peak shift)
#                 # , however the duration of the filter response is determined by the
#                 # max. period (e.g. 0.01-0.1Hz -> ~50s). So the time window chosen
#                 # should not be affected by the relatively small peak shift.
#                 # -- F * d
#                 obs_filt = signal.filtfilt(filter_b, filter_a, obs)
#                 # -- F * u (u = S*grn)
#                 syn_filt = signal.filtfilt(filter_b, filter_a, syn)
#
#                 # ------ apply window taper and polarity projection
#                 # obs = w * F * d
#                 obs_filt_win = np.dot(proj_matrix, obs_filt) * win_func
#                 # syn = w * F * u (u = S*g)
#                 syn_filt_win = np.dot(proj_matrix, syn_filt) * win_func
#
#                 # ------ measure CC time shift (between w*F*d and w*F*u)
#                 obs_norm = np.sqrt(np.sum(obs_filt_win**2))
#                 syn_norm = np.sqrt(np.sum(syn_filt_win**2))
#                 # window normalization factor (without dt)
#                 Nw = obs_norm * syn_norm
#                 # -- zero-lag normalized cc coeff.
#                 CC0 = np.sum(obs_filt_win * syn_filt_win)
#                 CC0 /= Nw
#
#                 # ------ measure adjoint source
#                 Aw = CC0 * obs_norm / syn_norm  # window amplitude raito
#                 if misfit_type == "cc0":
#                     # misfit: zero-lag cross-correlation
#                     # adjoint source: dchiw_du (misfit functional: zero-lag cc coef.)
#                     # dchiw_du = conj(F * [S]) * w * [ w * F * d - A * w * F * S * g] / N,
#                     # , where A = CC0(un-normalized) / norm(u)**2, N = norm(d)*norm(u)
#                     # -- dchiw_du
#                     # NOTE: *dt is put back to Nw
#                     dchiw_du1 = (
#                         win_func * (obs_filt_win - Aw * syn_filt_win) / Nw / syn_delta
#                     )
#                     # apply conj(F), equivalent to conj(F*conj(adj))
#                     # for two-pass filter (zero phase) conj(F) = F
#                     dchiw_du = signal.filtfilt(filter_b, filter_a, dchiw_du1[:, ::-1])
#                     dchiw_du = dchiw_du[:, ::-1]
#                 else:
#                     error_str = "%s:%s: unknown misfit type (%s)" % (
#                         station_id,
#                         window_id,
#                         misfit_type,
#                     )
#                     raise Exception(error_str)
#
#                 # add into total dchi_du
#                 dchi_du += weight * dchiw_du
#             # ------ end for window_id in windows:
#
#             # ====== output adjoint source for this station
#             # loop ENZ
#             for i in range(3):
#                 tr.data = dchi_du[i, syn_nl : (syn_nl + npts)]
#                 tr.stats.starttime = starttime
#                 tr.stats.delta = syn_delta
#
#                 out_file = "{:s}/{:s}.{:2s}{:1s}".format(
#                     out_dir, station_id, syn_band_code, syn_orientation_codes[i]
#                 )
#
#                 # sac format
#                 tr.write(out_file + ".adj.sac", "sac")
#
#                 # ascii format (needed by SEM)
#                 # time is relative to event origin time: t0
#                 with open(out_file + ".adj", "w") as fp:
#                     for j in range(npts):
#                         fp.write(
#                             "{:16.9e}  {:16.9e}\n".format(
#                                 syn_times[j], dchi_du[i, syn_nl + j]
#                             )
#                         )
#
#         # endfor station_id in station_dict:
#
#     # enddef measure_windows_for_one_station(self,
#
#
# ##
# ##======================================================
# ##
# #
# #  def output_adj_hess_part1(self,
# #      out_dir='adj',
# #      syn_band_code='MX'):
# #    """
# #    Output adjoint sources for the estimation of approximated Hessian diagonals for CC0 misfit.
# #
# #    This only output adjoint source of Part 1, that is, only the random part.
# #
# #    We need run adjoint simulations for Part 1 and 2 separately.
# #
# #    There are actually two more terms in the Hessian related to the recorded waveforms.
# #    However if the recorded and modelled waveforms are similar (cc0 close to 1), then
# #    the difference can be ignored (for hessian) and simplified to Part 2.
# #
# #    Notes
# #    -----
# #    For the objective function as the normalized zero-lag correlation ceofficient,
# #    the approximated Hessian can be seperated into two parts:
# #
# #    Part 1:
# #        H = - CC0 * norm(wFu)^(-2) * (wFdu1, wFdu2)
# #          = - Nw^-1 * Aw * (wFdu1, wFdu2)
# #        , where CC0 > 0 is the zero-lag correlation coefficient of (wFd, wFu)
# #        and the corresponding adjoint source is
# #
# #        adj = sqrt(CC0) * norm(wFu)^(-1) * conj(F)(transpose(w)r)
# #        , where r is a radome time series such that E(r*transpose(r)) = i.d.
# #        , for example, standard normal distribution N(0,1)
# #
# #        This follows the idea of random phase encoding method.
# #
# #    Part 2:
# #        H(x,y) = + CC0 * norm(wFu)^(-4) * (wFu, wFdu1) * (wFu, wFdu2)
# #               = H1(x) * H2(y)
# #        H1 = sqrt(CC0) * norm(wFu)^(-2) * (wFu, wFdu1)
# #        , the cooresponding adjoint source is
# #
# #        adj = sqrt(CC0) * norm(wFu)^(-2) * conj(F)(transpose(w)wFu)
# #    """
# #    syn_orientation_codes = ['E', 'N', 'Z']
# #    event = self.data['event']
# #
# #    #------ loop each station
# #    station_dict = self.data['station']
# #    for station_id in station_dict:
# #      station = station_dict[station_id]
# #      # skip rejected statations
# #      if station['stat']['code'] < 1:
# #        continue
# #
# #      # waveform
# #      waveform = station['waveform']
# #      time_sample = waveform['time_sample']
# #      syn_starttime = time_sample['starttime']
# #      syn_delta = time_sample['delta']
# #      syn_nt = time_sample['nt']
# #      syn_nl = time_sample['nl']
# #      syn_nr = time_sample['nr']
# #
# #      syn = waveform['syn']
# #
# #      #------ loop each window
# #      adj = np.zeros((3, syn_nt))
# #      # normal distribution
# #      rand = np.random.randn(3, syn_nt)
# #
# #      window_dict = station['window']
# #      for window_id in window_dict:
# #        # window parameters
# #        window = window_dict[window_id]
# #        # skip bad windows
# #        if window['stat']['code'] < 1:
# #          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
# #          continue
# #        if window['weight'] < 1.0e-3:
# #          continue
# #
# #        #------ window parameters
# #        # filter
# #        filter_dict = window['filter']
# #        filter_a = filter_dict['a']
# #        filter_b = filter_dict['b']
# #        # taper
# #        win_func = window['taper']['win']
# #        # polarity projection
# #        proj_matrix = window['polarity']['proj_matrix']
# #        # CC0
# #        cc0 = window['cc']['CC0']
# #
# #        #------ Part 1: random adjoint source
# #        # adj = sqrt(CC0) * norm(wFu)^(-1) * conj(F)(transpose(w)r)
# #        wr = np.dot(np.transpose(proj_matrix), rand) * win_func
# #        Fwr = signal.filtfilt(filter_b, filter_a, wr[:,::-1])
# #        Fwr = Fwr[:,::-1]
# #
# #        #------ Part 2: synthetic related adjoint source
# #        # adj = sqrt(CC0) * norm(wFu)^(-2) * conj(F)(transpose(w)wFu)
# #        Fu = signal.filtfilt(filter_b, filter_a, syn)
# #        wFu = np.dot(proj_matrix, Fu) * win_func
# #        norm_wFu = np.sqrt(np.sum(wFu**2))
# #        #wwFu = np.dot(np.transpose(proj_matrix), wFu) * win_func
# #        #FwwFu = signal.filtfilt(filter_b, filter_a, wwFu[:,::-1])
# #        #FwwFu = FwwFu[:,::-1]
# #
# #        #------ Part 1: make adjoint source for the current window
# #        adj_w = np.sqrt(cc0)/norm_wFu * Fwr
# #
# #        #------ Part 1: add into total adjoint source
# #        # chi = sum(chi_w * window_weight, over all w[indow])
# #        adj += np.sqrt(window['weight']) * adj_w
# #
# #      #endfor window_id in window_dict:
# #
# #      #------ output adjoint source
# #      # without padding
# #      npts = syn_nt - syn_nl - syn_nr
# #      starttime = syn_starttime + syn_nl*syn_delta
# #      # time samples for ascii output, referred to origin time
# #      syn_times = np.arange(npts)*syn_delta
# #      b = starttime - event['t0']
# #      syn_times += b
# #
# #      # loop ENZ
# #      tr = Trace()
# #      for i in range(3):
# #        tr.data = adj[i, syn_nl:(syn_nl+npts)]
# #        tr.stats.starttime = starttime
# #        tr.stats.delta = syn_delta
# #
# #        out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
# #            out_dir, station_id, syn_band_code,
# #            syn_orientation_codes[i])
# #
# #        # sac format
# #        tr.write(out_file + '.adj.sac', 'sac')
# #
# #        # ascii format (needed by SEM)
# #        # time is relative to event origin time: t0
# #        with open(out_file+'.adj','w') as fp:
# #          for j in range(npts):
# #            fp.write("{:16.9e}  {:16.9e}\n".format(
# #              syn_times[j], adj[i,syn_nl+j]))
# #
# #    #endfor station_id in station_dict:
# #
# ##
# ##======================================================
# ##
# #
# #  def output_adj_hess_part2(self,
# #      out_dir='adj',
# #      syn_band_code='MX'):
# #    """
# #    Output adjoint sources for the estimation of approximated Hessian diagonals for CC0 misfit.
# #
# #    This only output adjoint source of Part 2, which is related to the synthetics.
# #
# #    We need run adjoint simulations for Part 1 and 2 separately.
# #
# #    There are actually two more terms in the Hessian related to the recorded waveforms.
# #    However if the recorded and modelled waveforms are similar (cc0 close to 1), then
# #    the difference can be ignored (for hessian) and simplified to Part 2.
# #
# #    Notes
# #    -----
# #    For the objective function as the normalized zero-lag correlation ceofficient,
# #    the approximated Hessian can be seperated into two parts:
# #
# #    Part 1:
# #        H = - CC0 * norm(wFu)^(-2) * (wFdu1, wFdu2)
# #          = - Nw^-1 * Aw * (wFdu1, wFdu2)
# #        , where CC0 > 0 is the zero-lag correlation coefficient of (wFd, wFu)
# #        and the corresponding adjoint source is
# #
# #        adj = sqrt(CC0) * norm(wFu)^(-1) * conj(F)(transpose(w)r)
# #        , where r is a radome time series such that E(r*transpose(r)) = i.d.
# #        , for example, standard normal distribution N(0,1)
# #
# #        This follows the idea of random phase encoding method.
# #
# #    Part 2:
# #        H(x,y) = + CC0 * norm(wFu)^(-4) * (wFu, wFdu1) * (wFu, wFdu2)
# #               = H1(x) * H2(y)
# #        H1 = sqrt(CC0) * norm(wFu)^(-2) * (wFu, wFdu1)
# #        , the cooresponding adjoint source is
# #
# #        adj = sqrt(CC0) * norm(wFu)^(-2) * conj(F)(transpose(w)wFu)
# #    """
# #    syn_orientation_codes = ['E', 'N', 'Z']
# #    event = self.data['event']
# #
# #    #------ loop each station
# #    station_dict = self.data['station']
# #    for station_id in station_dict:
# #      station = station_dict[station_id]
# #      # skip rejected statations
# #      if station['stat']['code'] < 1:
# #        continue
# #
# #      # waveform
# #      waveform = station['waveform']
# #      time_sample = waveform['time_sample']
# #      syn_starttime = time_sample['starttime']
# #      syn_delta = time_sample['delta']
# #      syn_nt = time_sample['nt']
# #      syn_nl = time_sample['nl']
# #      syn_nr = time_sample['nr']
# #
# #      syn = waveform['syn']
# #
# #      #------ loop each window
# #      adj = np.zeros((3, syn_nt))
# #      # normal distribution
# #      #rand = np.random.randn(3, syn_nt)
# #
# #      window_dict = station['window']
# #      for window_id in window_dict:
# #        # window parameters
# #        window = window_dict[window_id]
# #        # skip bad windows
# #        if window['stat']['code'] < 1:
# #          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
# #          continue
# #        if window['weight'] < 1.0e-3:
# #          continue
# #
# #        #------ window parameters
# #        # filter
# #        filter_dict = window['filter']
# #        filter_a = filter_dict['a']
# #        filter_b = filter_dict['b']
# #        # taper
# #        win_func = window['taper']['win']
# #        # polarity projection
# #        proj_matrix = window['polarity']['proj_matrix']
# #        # CC0
# #        cc0 = window['cc']['CC0']
# #
# #        #------ Part 1: random adjoint source
# #        # adj = sqrt(CC0) * norm(wFu)^(-1) * conj(F)(transpose(w)r)
# #        #wr = np.dot(np.transpose(proj_matrix), rand) * win_func
# #        #Fwr = signal.filtfilt(filter_b, filter_a, wr[:,::-1])
# #        #Fwr = Fwr[:,::-1]
# #        #adj_w = np.sqrt(cc0)/norm_wFu * Fwr
# #
# #        #------ Part 2: synthetic related adjoint source
# #        # adj = sqrt(CC0) * norm(wFu)^(-2) * conj(F)(transpose(w)wFu)
# #        Fu = signal.filtfilt(filter_b, filter_a, syn)
# #        wFu = np.dot(proj_matrix, Fu) * win_func
# #        norm_wFu = np.sqrt(np.sum(wFu**2))
# #        wwFu = np.dot(np.transpose(proj_matrix), wFu) * win_func
# #        FwwFu = signal.filtfilt(filter_b, filter_a, wwFu[:,::-1])
# #        FwwFu = FwwFu[:,::-1]
# #        adj_w = np.sqrt(cc0) * norm_wFu**(-2) * FwwFu
# #
# #        #------ add into total adjoint source
# #        # chi = sum(chi_w * window_weight, over all w[indow])
# #        adj += np.sqrt(window['weight']) * adj_w
# #
# #      #endfor window_id in window_dict:
# #
# #      #------ output adjoint source
# #      # without padding
# #      npts = syn_nt - syn_nl - syn_nr
# #      starttime = syn_starttime + syn_nl*syn_delta
# #      # time samples for ascii output, referred to origin time
# #      syn_times = np.arange(npts)*syn_delta
# #      b = starttime - event['t0']
# #      syn_times += b
# #
# #      # loop ENZ
# #      tr = Trace()
# #      for i in range(3):
# #        tr.data = adj[i, syn_nl:(syn_nl+npts)]
# #        tr.stats.starttime = starttime
# #        tr.stats.delta = syn_delta
# #
# #        out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
# #            out_dir, station_id, syn_band_code,
# #            syn_orientation_codes[i])
# #
# #        # sac format
# #        tr.write(out_file + '.adj.sac', 'sac')
# #
# #        # ascii format (needed by SEM)
# #        # time is relative to event origin time: t0
# #        with open(out_file+'.adj','w') as fp:
# #          for j in range(npts):
# #            fp.write("{:16.9e}  {:16.9e}\n".format(
# #              syn_times[j], adj[i,syn_nl+j]))
# #
# #    #endfor station_id in station_dict:
# #
# #
# ##
# ##======================================================
# ##
# #
# #  def hess_diag_dmodel(self, model_name):
# #    """
# #    Output approximated Hessian diagonals for one model perturbation i.e. H(dm, dm).
# #
# #    There are actually two more terms in the Hessian related to the recorded waveforms.
# #    However if the recorded and modelled waveforms are similar (cc0 close to 1), then
# #    the difference can be ignored (for hessian) and simplified to Part 2.
# #
# #    Notes
# #    -----
# #    For the objective function as the normalized zero-lag correlation ceofficient,
# #    the approximated Hessian can be seperated into two parts:
# #
# #    Part 1:
# #        H = - CC0 * norm(wFu)^(-2) * (wFdu1, wFdu2)
# #          = - Nw^-1 * Aw * (wFdu1, wFdu2)
# #        , where CC0 > 0 is the zero-lag correlation coefficient of (wFd, wFu)
# #
# #    Part 2:
# #        H(x,y) = + CC0 * norm(wFu)^(-4) * (wFu, wFdu1) * (wFu, wFdu2)
# #               = H1(x) * H2(y)
# #    """
# #    event = self.data['event']
# #
# #    #------ loop each station
# #    station_dict = self.data['station']
# #    for station_id in station_dict:
# #      station = station_dict[station_id]
# #      # skip rejected statations
# #      if station['stat']['code'] < 1:
# #        continue
# #
# #      # initial waveform
# #      waveform = station['waveform']
# #      time_sample = waveform['time_sample']
# #      syn_starttime = time_sample['starttime']
# #      syn_delta = time_sample['delta']
# #      syn_nt = time_sample['nt']
# #      syn_nl = time_sample['nl']
# #      syn_nr = time_sample['nr']
# #      syn = waveform['syn']
# #
# #      # perturbed waveform
# #      waveform_der = station['waveform_der']
# #      du = waveform_der[model_name]['du']
# #
# #      #------ loop each window
# #      window_dict = station['window']
# #      for window_id in window_dict:
# #        # window parameters
# #        window = window_dict[window_id]
# #        # skip bad windows
# #        if window['stat']['code'] < 1:
# #          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
# #          continue
# #
# #        #------ window parameters
# #        # filter
# #        filter_dict = window['filter']
# #        filter_a = filter_dict['a']
# #        filter_b = filter_dict['b']
# #        # taper
# #        win_func = window['taper']['win']
# #        # polarity projection
# #        proj_matrix = window['polarity']['proj_matrix']
# #        # CC0
# #        cc0 = window['cc']['CC0']
# #
# #        #------ hessian for current window
# #        Fu = signal.filtfilt(filter_b, filter_a, syn)
# #        wFu = np.dot(proj_matrix, Fu) * win_func
# #        norm2_wFu = np.sum(wFu**2)
# #
# #        Fdu = signal.filtfilt(filter_b, filter_a, du)
# #        wFdu = np.dot(proj_matrix, Fdu) * win_func
# #
# #        hess_win = -1.0 * cc0/norm2_wFu * (np.sum(wFdu**2) - np.sum(wFu*wFdu)**2/norm2_wFu)
# #
# #        #------ record hess
# #        if 'hess_diag' not in window:
# #          window['hess_diag'] = {}
# #        window['hess_diag'][model_name] = hess_win
# #
# #      #endfor window_id in window_dict:
# #    #endfor station_id in station_dict:
# #  #enddef hess_diag_dmodel
# #
# ##
# ##======================================================
# ##
# #
# #  def output_hess_diag(self, model_name, out_file='hess_diag.txt'):
# #    """
# #    Output hess diag
# #
# #    """
# #    event = self.data['event']
# #    station_dict = self.data['station']
# #
# #    f = open(out_file, 'w')
# #    f.write("#station window hess_diag(%s)\n" % (model_name))
# #
# #    #------ loop each station
# #    for station_id in station_dict:
# #      station = station_dict[station_id]
# #      # skip rejected statations
# #      if station['stat']['code'] < 1:
# #        continue
# #
# #      #------ loop each window
# #      window_dict = station['window']
# #      for window_id in window_dict:
# #        # window parameters
# #        window = window_dict[window_id]
# #        # skip bad windows
# #        if window['stat']['code'] < 1:
# #          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
# #          continue
# #        # write out hess_diag
# #        hess_diag = window['hess_diag'][model_name]
# #
# #        f.write("{:15s} {:15s} {:12.5e}\n".format(
# #          station_id, window_id, hess_diag))
# #
# #    f.close()
# #  #enddef
# #
# ##
# ##======================================================
# ##
# #
# #  def output_adj_hess_model_product(self,
# #      model_name='perturb',
# #      out_dir='adj',
# #      syn_band_code='MX'):
# #    """
# #    Output adjoint sources for calculating Hessian-model product H * dm.
# #
# #    Notes
# #    -----
# #    For the objective function as the normalized zero-lag correlation ceofficient,
# #    the approximate Hessian for one data window is (ignore second order derivative in u)
# #
# #    Hw = - Nw^-1 * Aw * (wFdu1, wFdu2)
# #         + 3 * Nw^-1 * Aw * norm(wFu)^-2 * (wFu, wFdu1) * (wFu, wFdu2)
# #         - Nw^-1 * norm(wFu)^-2 * (wFu, wFdu1) * (wFd, wFdu2)
# #         - Nw^-1 * norm(wFu)^-2 * (wFu, wFdu2) * (wFd, wFdu1)
# #
# #    , where Nw = norm(wFu)*norm(wFd), Aw = (wFu, wFd)/norm(wFu)^2
# #    and norm(.) = sqrt((., .)), (.,.) is inner product. Notice Aw/Nw = cc0/norm(wFu)^2
# #
# #    If the recorded and modelled waveforms are close enough (cc0 close to 1), then
# #    the last three terms in Hw can be simplified to one term (wFd ~ Aw * wFu):
# #
# #        + Nw^-1 * Aw * norm(wFu)^-2 * (wFu, wFdu1) * (wFu, wFdu2)
# #
# #    For the Hessian-model product H(., dm) let du2 = Du(m; dm) = u(m+dm) - u(m)
# #    then the adjoint source is (simplified)
# #
# #    r = - Nw^-1 * Aw * wFdu2 + Nw^-1 * Aw * norm(wFu)^-2 * (wFu, wFdu2) * wFu
# #
# #    adj_w = conj(F)(transpose(w) * r)
# #
# #    The total hessian is sum(weight * Hw) and the total adjoint source is sum(weight * adj_w)
# #
# #    """
# #    syn_orientation_codes = ['E', 'N', 'Z']
# #    event = self.data['event']
# #
# #    #------ loop each station
# #    station_dict = self.data['station']
# #    for station_id in station_dict:
# #      station = station_dict[station_id]
# #      # skip rejected statations
# #      if station['stat']['code'] < 1:
# #        continue
# #
# #      # waveform
# #      waveform = station['waveform']
# #      time_sample = waveform['time_sample']
# #      syn_starttime = time_sample['starttime']
# #      syn_delta = time_sample['delta']
# #      syn_nt = time_sample['nt']
# #      syn_nl = time_sample['nl']
# #      syn_nr = time_sample['nr']
# #
# #      syn = waveform['syn']
# #
# #      du = station['waveform_der'][model_name]['du']
# #
# #      #------ loop each window
# #      adj = np.zeros((3, syn_nt))
# #
# #      window_dict = station['window']
# #      for window_id in window_dict:
# #        # window parameters
# #        window = window_dict[window_id]
# #        # skip bad windows
# #        if window['stat']['code'] < 1:
# #          warnings.warn("Window %s not measured for adj, SKIP" % window_id)
# #          continue
# #        if window['weight'] < 1.0e-3:
# #          continue
# #
# #        #------ window parameters
# #        # filter
# #        filter_dict = window['filter']
# #        filter_a = filter_dict['a']
# #        filter_b = filter_dict['b']
# #        # taper
# #        win_func = window['taper']['win']
# #        # polarity projection
# #        proj_matrix = window['polarity']['proj_matrix']
# #        # CC0
# #        cc0 = window['cc']['CC0']
# #
# #        #------ make adjoint source
# #        # r = - Nw^-1 * Aw * wFdu + Nw^-1 * Aw * norm(wFu)^-2 * (wFu, wFdu) * wFu
# #        #   = - cc0/norm(wFu)^2 * (wFdu - norm(wFu)^-2*(wFu, wFdu)*wFu)
# #        # adj_w = conj(F)(transpose(w) * r)
# #        Fu = signal.filtfilt(filter_b, filter_a, syn)
# #        wFu = np.dot(proj_matrix, Fu) * win_func
# #        norm_wFu = np.sqrt(np.sum(wFu**2))
# #
# #        Fdu = signal.filtfilt(filter_b, filter_a, du)
# #        wFdu = np.dot(proj_matrix, Fdu) * win_func
# #
# #        adj_w = -cc0/norm_wFu**2 * (wFdu - wFu*np.sum(wFu*wFdu)/norm_wFu**2)
# #        adj_w = np.dot(np.transpose(proj_matrix), adj_w) * win_func
# #        adj_w = signal.filtfilt(filter_b, filter_a, adj_w[:,::-1])
# #        adj_w = adj_w[:,::-1]
# #
# #        #------ add into total adjoint source
# #        adj += window['weight'] * adj_w
# #
# #      #endfor window_id in window_dict:
# #
# #      #------ output adjoint source
# #      # without padding
# #      npts = syn_nt - syn_nl - syn_nr
# #      starttime = syn_starttime + syn_nl*syn_delta
# #      # time samples for ascii output, referred to origin time
# #      syn_times = np.arange(npts)*syn_delta
# #      b = starttime - event['t0']
# #      syn_times += b
# #
# #      # loop ENZ
# #      tr = Trace()
# #      for i in range(3):
# #        tr.data = adj[i, syn_nl:(syn_nl+npts)]
# #        tr.stats.starttime = starttime
# #        tr.stats.delta = syn_delta
# #
# #        out_file = '{:s}/{:s}.{:2s}{:1s}'.format(
# #            out_dir, station_id, syn_band_code,
# #            syn_orientation_codes[i])
# #
# #        # sac format
# #        tr.write(out_file + '.adj.sac', 'sac')
# #
# #        # ascii format (needed by SEM)
# #        # time is relative to event origin time: t0
# #        with open(out_file+'.adj','w') as fp:
# #          for j in range(npts):
# #            fp.write("{:16.9e}  {:16.9e}\n".format(
# #              syn_times[j], adj[i,syn_nl+j]))
# #
# #    #endfor station_id in station_dict:
#
# # END class misfit
