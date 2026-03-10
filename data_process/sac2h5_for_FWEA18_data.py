#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert SEM sac files to h5 format with filtering and resampling.

This module provides functionality to convert seismic waveform data from SAC format
to HDF5 format, applying low-pass filtering and resampling as needed.
"""

import os
import argparse
import warnings
from typing import List, Optional

import scipy
import scipy.signal
import scipy.fft
import numpy as np
import tables
import pandas as pd
import obspy
from obspy import Stream, Trace, read_events


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="convert SEM sac files to h5")

    parser.add_argument("--sac_dir", required=True, help="directory of SAC files")
    parser.add_argument("--channel_file", required=True, help="channel file")
    parser.add_argument("--out_h5", required=True, help="output h5 file")
    parser.add_argument("--cmt_file", default=None, help="CMTSOLUTION file")
    parser.add_argument(
        "--resample_fs", type=float, default=2.0, help="resampling frequency"
    )
    parser.add_argument(
        "--low_pass_freq",
        type=float,
        default=0.5,
        help="corner frequency of butterworth low pass filter",
    )
    parser.add_argument(
        "--low_pass_order",
        type=int,
        default=15,
        help="order of the butterworth low pass filter",
    )
    parser.add_argument(
        "--high_pass_freq",
        type=float,
        default=0.01,
        help="corner frequency of butterworth high pass filter",
    )
    parser.add_argument(
        "--high_pass_order",
        type=int,
        default=5,
        help="order of the butterworth high pass filter",
    )
    parser.add_argument("--sac_suffix", default="", help="suffix of SEM sac files")
    parser.add_argument(
        "--data_type", default="DISP", help="data type of SEM sac files, e.g. VEL, DISP"
    )

    return parser.parse_args()


def apply_filter(
    trace: Trace,
    low_pass_freq: float,
    low_pass_order: int,
    high_pass_freq: float,
    high_pass_order: int,
) -> None:
    """
    Apply low-pass filter to a trace using a Butterworth filter.

    Args:
        trace: ObsPy trace to filter
        low_pass_freq: Corner frequency for the low-pass filter
        low_pass_order: Order of the Butterworth filter
    """
    fs = trace.stats.sampling_rate
    npts = trace.stats.npts

    # Calculate padding and FFT length
    pad_length = max(10.0, 2.0 / low_pass_freq)
    npad = int(pad_length * fs)
    nfft = scipy.fft.next_fast_len(npts + 2 * npad)
    freqs = np.fft.rfftfreq(nfft, d=1 / fs)

    # create low-pass filter in frequency domain
    lp_sos = scipy.signal.butter(
        low_pass_order, low_pass_freq, "lowpass", fs=fs, output="sos"
    )
    _, h_lp = scipy.signal.freqz_sos(lp_sos, worN=freqs, fs=fs)

    # create high-pass filter in frequency domain
    hp_sos = scipy.signal.butter(
        high_pass_order, high_pass_freq, "highpass", fs=fs, output="sos"
    )
    _, h_hp = scipy.signal.freqz_sos(hp_sos, worN=freqs, fs=fs)

    # Transform to frequency domain, apply filter, transform back
    sig_spectrum = np.fft.rfft(
        np.pad(trace.data, (npad, nfft - npts - npad), mode="edge"), nfft
    )
    sig_spectrum *= abs(h_lp) * abs(h_hp)
    filtered_data = np.fft.irfft(sig_spectrum, nfft)[npad : (npts + npad)]

    trace.data = filtered_data


def interpolate_trace(
    trace: Trace, resample_fs: float, resample_starttime: float, resample_npts: int
) -> None:
    """
    Interpolate trace to desired sampling rate and time window.

    Args:
        trace: ObsPy trace to interpolate
        resample_fs: Target sampling frequency
        resample_starttime: Start time for resampling
        resample_npts: Number of points in resampled trace
    """
    trace.interpolate(
        resample_fs,
        starttime=resample_starttime,
        npts=resample_npts,
        method="lanczos",
        a=20,
    )


def load_channel_data(channel_file: str) -> pd.DataFrame:
    """
    Load station information from CSV file.

    Args:
        station_file: Path to station file

    Returns:
        DataFrame with station information
    """
    df = pd.read_csv(
        channel_file,
        sep=r"|",
        header=None,
        dtype=str,
        keep_default_na=False,
        usecols=range(10),
    )
    df.columns = [
        "net",
        "sta",
        "loc",
        "cha",
        "lat",
        "lon",
        "elevation",
        "depth",
        "azimuth",
        "dip",
    ]
    return df


def process_station_traces(
    sac_dir: str,
    net: str,
    sta: str,
    loc: str,
    channels: List[str],
    sac_suffix: str,
    resample_fs: float,
    low_pass_freq: float,
    low_pass_order: int,
    high_pass_freq: float,
    high_pass_order: int,
) -> Optional[Stream]:
    """
    Process all traces for a single station.

    Args:
        sac_dir: Directory containing SAC files
        net: Network name
        sta: Station name
        loc: Location code
        channels: List of channel names to process
        sac_suffix: Suffix for SAC files
        resample_fs: Resampling frequency
        low_pass_freq: Low-pass filter frequency
        low_pass_order: Low-pass filter order
        pad_before: time to pad before each trace (seconds)

    Returns:
        ObsPy Stream object with processed traces, or None if failed
    """
    st = Stream()

    # Load all channel traces
    for cha in channels:
        sac_file = os.path.join(sac_dir, f"{net}.{sta}.{loc}.{cha}{sac_suffix}")
        if not os.path.exists(sac_file):
            warnings.warn(
                f"[WARN] {sac_file} does not exist, skipping station {net}.{sta}."
            )
            return None

        try:
            tr = obspy.read(sac_file)
            st.extend(tr)
        except Exception as e:
            warnings.warn(
                f"[WARN] Failed to read {sac_file}: {e}, skipping station {net}.{sta}."
            )
            return None

    # Determine common time window for resampling
    max_starttime = min(tr.stats.starttime for tr in st)
    min_endtime = max(tr.stats.endtime for tr in st)
    resample_starttime = max_starttime
    resample_npts = int((min_endtime - resample_starttime) * resample_fs)
    # pad_npts, pad_time = 0, 0.0
    # if pad_before > 0:
    #     pad_npts = int(pad_before * resample_fs)
    #     pad_time = pad_npts / resample_fs
    resample_endtime = resample_starttime + (resample_npts - 1) / resample_fs

    # Process each trace: filter and resample
    for tr in st:
        npad = int((tr.stats.starttime - resample_starttime) * tr.stats.sampling_rate)
        # add extra padding to ensure interpolation within the old time range
        npad_before = max(5, npad)

        npad = int((resample_endtime - tr.stats.endtime) * tr.stats.sampling_rate)
        npad_after = max(5, npad)

        tr.data = np.pad(tr.data, (npad_before, npad_after), mode="edge")
        tr.stats.starttime = tr.stats.starttime - npad_before / tr.stats.sampling_rate

        apply_filter(tr, low_pass_freq, low_pass_order, high_pass_freq, high_pass_order)
        interpolate_trace(tr, resample_fs, resample_starttime, resample_npts)

    return st


def create_h5_dataset(
    h5f: tables.File,
    group_name: str,
    array_name: str,
    data: np.ndarray,
    attrs: dict,
    filters: tables.Filters,
) -> tables.CArray:
    """
    Create HDF5 dataset with specified attributes.

    Args:
        h5f: HDF5 file handle
        group_name: Name of the group to create dataset in
        data: Data array to store
        attrs: Attributes dictionary
        filters: HDF5 compression filters

    Returns:
        Created CArray object
    """
    atom = tables.Atom.from_dtype(np.dtype(np.float32))
    # array_name = attrs.get("array_name", "DATA")

    # Create dataset
    ca = h5f.create_carray(
        h5f.root[group_name], array_name, atom, data.shape, filters=filters
    )
    ca[:] = data.astype(np.float32)

    # Set attributes
    for key, value in attrs.items():
        setattr(ca.attrs, key, value)

    return ca


def sac2h5(
    sac_dir: str,
    channel_file: str,
    h5_file: str,
    cmt_file: Optional[str] = None,
    resample_fs: float = 2.0,
    low_pass_freq: float = 0.5,
    low_pass_order: int = 15,
    high_pass_freq: float = 0.01,
    high_pass_order: int = 5,
    sac_suffix: str = ".sem.sac",
    data_type: str = "DISP",
):
    """
    Convert SAC files to HDF5 format with filtering and resampling.

    Args:
        sac_dir: Directory containing SAC files
        station_file: File with station information
        cmt_file: CMTSOLUTION file
        resample_fs: Target sampling frequency
        low_pass_freq: Low-pass filter corner frequency
        low_pass_order: Low-pass filter order
        channels: List of channel names to process
        sac_suffix: Suffix for SAC files
        data_type: Data type (DISP, VEL, etc.)
        h5_file: Output HDF5 file path
    """
    # Open HDF5 file
    with tables.open_file(h5_file, mode="w") as h5f:
        groot = h5f.root

        if cmt_file is not None:
            event_data = read_events(cmt_file)[0]
            groot._v_attrs["event"] = event_data

        channel_data = load_channel_data(channel_file)

        for key, channels in channel_data.groupby(["net", "sta", "loc"]):
            net, sta, loc = key

            print(f"[INFO] processing {net}.{sta}.{loc}")

            lats = channels["lat"].astype(float).unique()
            if len(lats) > 1:
                raise ValueError(
                    f"[ERROR] multiple latitudes for {net}.{sta}.{loc}: {lats}"
                )
            stla = lats[0]

            lons = channels["lon"].astype(float).unique()
            if len(lons) > 1:
                raise ValueError(
                    f"[ERROR] multiple longitudes for {net}.{sta}.{loc}: {lons}"
                )
            stlo = lons[0]

            elevs = channels["elevation"].astype(float).unique()
            if len(elevs) > 1:
                raise ValueError(
                    f"[ERROR] multiple elevations for {net}.{sta}.{loc}: {elevs}"
                )
            stel = elevs[0]

            depths = channels["depth"].astype(float).unique()
            if len(depths) > 1:
                raise ValueError(
                    f"[ERROR] multiple depths for {net}.{sta}.{loc}: {depths}"
                )
            stdp = depths[0]

            # Process all traces for this station
            st = process_station_traces(
                sac_dir,
                net,
                sta,
                loc,
                list(channels["cha"]),
                sac_suffix,
                resample_fs,
                low_pass_freq,
                low_pass_order,
                high_pass_freq,
                high_pass_order,
            )

            if st is None:
                continue  # Skip if traces couldn't be processed

            assert len(st) == len(channels)

            # Create or overwrite station group
            station_name = f"{net}_{sta}"
            if station_name in groot:
                warnings.warn(f"[WARN] {station_name} already exists, overwriting.")
                h5f.remove_node(groot, station_name)

            gsta = h5f.create_group(groot, station_name)
            gsta._v_attrs["network"] = net
            gsta._v_attrs["station"] = sta
            gsta._v_attrs["location"] = loc
            gsta._v_attrs["latitude"] = float(stla)
            gsta._v_attrs["longitude"] = float(stlo)
            gsta._v_attrs["elevation"] = float(stel)
            gsta._v_attrs["depth"] = float(stdp)

            # Prepare waveform data array
            nchan = len(st)
            npts = st[0].stats.npts
            # resample_npts = int((max(tr.stats.endtime for tr in st) -
            #                     min(tr.stats.starttime for tr in st)) * resample_fs)
            shape = (nchan, npts)

            # Create compressed dataset
            compression_filter = tables.Filters(complevel=3, complib="zlib")
            waveform_data = np.zeros(shape, dtype=np.float32)
            for i in range(nchan):
                waveform_data[i, :] = np.array(st[i].data, dtype=np.float32)

            # Define channel metadata
            dtype = [("name", "S3"), ("azimuth", float), ("dip", float)]
            channel_info = [
                (
                    cha["cha"],
                    float(cha["azimuth"]),
                    float(cha["dip"]),
                )
                for _, cha in channels.iterrows()
            ]

            # Create dataset and set attributes
            array_name = f"DATA_{data_type}"
            attrs = {
                "channels": np.array(channel_info, dtype=dtype),
                "starttime": st[0].stats.starttime,
                "sampling_rate": resample_fs,
                "npts": npts,
                "filter": [
                    {"type": "lowpass", "N": low_pass_order, "Wn": low_pass_freq},
                    {"type": "highpass", "N": high_pass_order, "Wn": high_pass_freq},
                ],
                "response_corner_frequency": [],
                "type": data_type,
            }

            create_h5_dataset(
                h5f, station_name, array_name, waveform_data, attrs, compression_filter
            )


def main():
    """Main entry point."""
    args = parse_arguments()
    print(f"Arguments: {args}")

    assert args.resample_fs > 0
    assert args.low_pass_freq > 0
    assert args.low_pass_order > 0
    assert args.low_pass_freq < args.resample_fs / 2
    assert args.high_pass_freq > 0
    assert args.high_pass_order > 0
    assert args.high_pass_freq < args.low_pass_freq

    sac2h5(
        args.sac_dir,
        args.channel_file,
        args.out_h5,
        cmt_file=args.cmt_file,
        resample_fs=args.resample_fs,
        low_pass_freq=args.low_pass_freq,
        low_pass_order=args.low_pass_order,
        high_pass_freq=args.high_pass_freq,
        high_pass_order=args.high_pass_order,
        sac_suffix=args.sac_suffix,
        data_type=args.data_type,
    )


if __name__ == "__main__":
    main()
