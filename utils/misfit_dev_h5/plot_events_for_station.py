#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""get model update step length from grid search results"""
import warnings
import argparse

import tables as pt
import numpy as np
import scipy

import pyproj
# from obspy.geodetics import gps2dist_azimuth
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.imaging.beachball import beach

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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


def _extract_obs_syn_ENZ(
    obs_data, obs_attrs, syn_data, syn_attrs  # , event_tau
):
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
        msg = f"problematic Z channel info: {obs_Zchan}, skip"
        raise AssertionError(msg)
    no_Hchan = False
    if len(obs_Hchan) == 0:
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

    return obs_ENZ, syn_ENZ, no_Hchan


def plot_seismogram_1comp(
    misfit_h5file,
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

    h5f = pt.open_file(misfit_h5file, "r")

    network = h5f.root._v_attrs["network"]
    station = h5f.root._v_attrs["station"] 
    stnm = f"{network}_{station}"

    # config
    config = h5f.root._v_attrs["config"]
    geod = pyproj.Geod(ellps=config["gps_ellps"])
    R_earth = (geod.a + geod.b) / 2
    taup_model = TauPyModel(model=config["taup_model"])
    obs_tag = config["data"]["tag"]
    syn_tag = config["syn"]["tag"]

    if "/waveform" not in h5f:
        msg = '"/waveform" does not exist!'
        raise KeyError(msg)
    g_wav = h5f.get_node("/waveform")

    if "/window" not in h5f:
        msg = '"/waveform" does not exist!'
        raise KeyError(msg)
    tbl_win = h5f.get_node("/window")

    # get windows
    windows = []
    for g_evt in g_wav:
        evnm = g_evt._v_attrs["event_name"]
        win = tbl_win.read_where(f'(id == b"{win_id}") & (event_name == b"{evnm}")')
        if len(win) == 0:
            continue
        win = win[0]
        # events.append((g_evt, win[0]))
        evla = g_evt._v_attrs["evla"]
        evlo = g_evt._v_attrs["evlo"]
        evdp = g_evt._v_attrs["evdp"]
        stla = g_evt._v_attrs["stla"]
        stlo = g_evt._v_attrs["stlo"]
        mt = g_evt._v_attrs["moment_tensor"]
        Mrr = mt[0][0]
        Mtt = mt[1][1]
        Mpp = mt[2][2]
        Mrt = mt[0][1]
        Mrp = mt[0][2]
        Mtp = mt[1][2]
        focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
        evt0 = UTCDateTime(g_evt._v_attrs["origin_time"])
        evtau = g_evt._v_attrs["evtau"]
        winb = UTCDateTime(win["starttime"]) - evt0
        wine = UTCDateTime(win["endtime"]) - evt0
        winc = (winb + wine) / 2.0
        az, baz, dist = geod.inv(evlo, evla, stlo, stla)
        az = az % 360
        baz = baz % 360
        dist_degree = np.rad2deg(dist / R_earth)
        windows.append(
            { 
                "win":          win, 
                "evnm":         evnm,
                "evla":         evla,
                "evlo":         evlo,
                "evdp":         evdp,
                "evt0":         evt0,
                "evtau":        evtau,
                "focmec":      focmec,
                "stla":         stla,
                "stlo":         stlo,
                "winb":         winb,
                "wine":         wine,
                "winc":         winc,
                "az":           baz,
                # "baz":          baz,
                "dist_degree":  dist_degree,
                "g_evt":        g_evt
            }
        )
    stla_all = np.array([w["stla"] for w in windows])
    stlo_all = np.array([w["stlo"] for w in windows])
    dist_all = np.array([w["dist_degree"] for w in windows])
    winc_all = np.array([w["winc"] for w in windows])

    # get average moveout of the window center
    # linear regression tc = dist*rayp + tb
    A = np.vstack([dist_all, np.ones(len(dist_all))]).T
    plot_rayp, plot_c = np.linalg.lstsq(A, winc_all, rcond=None)[0]
    # round to integer
    plot_rayp = np.round(plot_rayp)

    # ------ calculate traveltime curves (only for body wave)
    phases = set([w["phase"].decode() for e in windows if (w := e["win"])["type"] == b"body"])
    phase_list = [a for p in phases for a in p.split(",")]
    if windows[0]["win"]["type"] == b"Pwin":
        phase_list.extend(["p", "P"])
    if windows[0]["win"]["type"] == b"Swin":
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
    windows = sorted(windows, key=lambda w: w["az"] % 360)  # sort windows by back-azimuth
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
            win = windows[bin_idx0]
            azmin = win["az"] % 360
            azmax = azmin + plot_azbin

            # get windows of one azimuthal bin
            windows_bin = []
            nwin_bin = 0
            while bin_idx0 < nwin:
                win = windows[bin_idx0]
                winfo = win["win"]
                if not winfo["valid"]:
                    msg = f"invalid window: {winfo}, skip"
                    warnings.warn(msg)
                    continue
                az = win["az"] % 360
                dist_degree = win["dist_degree"]
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
                # put window into current bin or start a new bin
                if az <= azmax and (
                    not plot_max_ntrace_per_bin
                    or nwin_bin <= plot_max_ntrace_per_bin
                ):
                    windows_bin.append(windows[bin_idx0])
                    nwin_bin = nwin_bin + 1
                    bin_idx0 += 1  # next window
                else:
                    break

            # skip empty azimuthal_bin
            if nwin_bin == 0:
                msg = f"No traces in the azimuthal bin [{bin_azmin}, {bin_azmax}], skip"
                warnings.warn(msg)
                continue

            bin_azmin = windows_bin[0]["az"] % 360
            bin_azmax = windows_bin[-1]["az"] % 360

            # ---- create figure
            print(
                f"[INFO] plot azimuthal bin: {bin_azmin:05.1f} - {bin_azmax:05.1f}"
            )

            # fig = plt.figure(figsize=(8.5, 11)) # US letter
            fig = plt.figure(figsize=(8.27, 11.69))  # A4
            str_title = "{:s}.{:s} ({:s} az:{:04.1f}~{:04.1f})".format(
                network, station, plot_window_id, bin_azmin, bin_azmax
            )
            if title_prefix:
                str_title = "{:s} {:s}".format(title_prefix, str_title)
            fig.suptitle(str_title)
            # fig.text(
            #    0.5, 0.965, str_title, size="x-large", horizontalalignment="center"
            # )

            stla_bin = np.array([win["stla"] for win in windows_bin])
            stlo_bin = np.array([win["stlo"] for win in windows_bin])
            # weight_bin = np.array([win["win"]["weight"] for win in windows_bin])
            # ccdt_bin = np.array([win["win"]["cc_time_shift"] for win in windows_bin])
            # cc0_bin = np.array([win["win"]["cc0"] for win in windows_bin])

            # event,station map with cc_dt
            ax_size = [0.35, 0.4]
            ax_origin = [0.02, 0.55]
            ax = fig.add_axes(ax_origin + ax_size, projection=projection)
            ax.set_extent(map_extent) #, crs=projection)
            # The facecolor argument controls the fill color
            # ax.add_feature(cfeature.OCEAN, facecolor='#AED6F1', zorder=0)
            # ax.add_feature(cfeature.LAND, facecolor='#FAD7A0', zorder=0)
            # Optional: Add coastlines on top for clarity
            ax.coastlines(resolution='50m', color='black', linewidth=0.5)
            ax.add_feature(cfeature.BORDERS, linewidth=0.2, alpha=0.5)
            # ax.set_xticks([10, 20, 30, 40, 50], crs=projection)
            # lat_formatter = LatitudeFormatter(number_format='.0f', 
            #                                   degree_symbol='', 
            #                                   direction_label=False)
            # lon_formatter = LongitudeFormatter(number_format='.0f', 
            #                                    degree_symbol='', 
            #                                    direction_label=False)
            ax.gridlines(draw_labels={"bottom": "x", "left": "y"}, 
                         # xformatter=lon_formatter, 
                         # yformatter=lat_formatter,
                         formatter_kwargs={'number_format': '.0f', 
                                           'degree_symbol': '', 
                                           'direction_label': False},
                         xlabel_style={'fontsize': 7},
                         ylabel_style={'fontsize': 7},
                         xpadding=1.5,
                         ypadding=1.5,
                         linewidth=0.2,)
            # ax.stock_img()
            # ax.set_title(f"lat:{evla:.3f} lon:{evlo:.3f} dep:{evdp:.1f}")

            # plot focal mechanism
            for win in windows:
                evla, evlo, focmec = win["evla"], win["evlo"], win["focmec"]
                sx, sy = projection.transform_point(evlo, evla, ccrs.Geodetic())
                x0, x1 = ax.get_xlim()
                bb_width = 0.04 * abs(x1 - x0)
                facecolor = "k"
                if win in windows_bin:
                    facecolor = "r"
                b = beach(
                    focmec, xy=(sx, sy), width=bb_width, linewidth=0.2, facecolor=facecolor, zorder=2
                )
                ax.add_collection(b)

            im = ax.scatter(
                stlo_bin,
                stla_bin,
                s=12,
                marker="^",
                facecolor="k",
                zorder=2,
                transform=ccrs.Geodetic(),
            )

            # create axis for seismograms
            ax_origin = [0.42, 0.06]
            ax_size = [0.4, 0.88]
            # ax_size = [0.3, 0.90]
            ax_1comp = fig.add_axes(ax_origin + ax_size)

            # -- xlim setting
            winb_bin = np.array([win["winb"] for win in windows_bin])
            wine_bin = np.array([win["wine"] for win in windows_bin])
            dist_bin = np.array([win["dist_degree"] for win in windows_bin])
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
            for data in windows_bin:
                win = data["win"]
                evnm = data["evnm"]
                evt0 = data["evt0"]
                evtau = data["evtau"]
                dist_degree = data["dist_degree"]

                g_sta = data["g_evt"]

                if obs_tag in g_sta:
                    obs = g_sta[obs_tag]
                elif "DATA_DISP" in g_sta: 
                    obs = g_sta["DATA_DISP"]
                elif "DATA_VEL" in g_sta: 
                    obs = g_sta["DATA_VEL"]
                else:
                    msg = f"{stnm}: {obs_tag}/DATA_DISP/DATA_VEL not in {g_sta}"
                    warnings.warn(msg)
                    continue

                if syn_tag not in g_sta:
                    msg = f"{stnm}: {syn_tag} not in {g_sta}"
                    warnings.warn(msg)
                    continue
                syn = g_sta[syn_tag]

                obs_attrs = {k: obs.attrs[k] for k in obs.attrs._f_list("user")}
                syn_attrs = {k: syn.attrs[k] for k in syn.attrs._f_list("user")}
                try:
                    obs_ENZ, syn_ENZ, no_Hchan = _extract_obs_syn_ENZ(
                        obs[:], obs_attrs, syn[:], syn_attrs
                    )
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
                str_annot = "%s,%.2f,%.1f,%.1f,%.1f" % (
                    evnm,
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
                    fontsize=5,
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


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("misfit_h5file")  # list of misfit.h5 of each events
    parser.add_argument("window_id")  # list of misfit.h5 of each events
    parser.add_argument("out_fig", help="output pdf file")
    parser.add_argument("--azbin_size", default=360, type=float, help="max szie of each azimuthal bin")
    parser.add_argument("--max_traces", default=100, type=int, help="max number of traces per azimuthal bin")
    # parser.add_argument("--parfile", default=None, help="re-read ymal file")
    args = parser.parse_args()

    plot_seismogram_1comp(args.misfit_h5file, 
                          args.window_id, 
                          out_fig=args.out_fig,
                          azbin=args.azbin_size, 
                          max_ntrace_per_bin=args.max_traces,)


if __name__ == "__main__":
    main()