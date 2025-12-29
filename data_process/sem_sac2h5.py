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
from typing import List, Optional, Tuple

import scipy
import scipy.signal
import scipy.fft
import numpy as np
import tables
import pandas as pd
import obspy
from obspy import Stream, Trace


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="convert SEM sac files to h5")

    parser.add_argument("sac_dir", help="directory of SAC files")
    parser.add_argument("station", help="SEM STATIONS file")
    parser.add_argument("out_h5", help="output h5 file")
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
        "--channels",
        nargs="+",
        default=["BXE", "BXN", "BXZ"],
        help="channel names of the SEM sac files, sac_dir/net.sta.[channel][sac_suffix]",
    )
    parser.add_argument(
        "--sac_suffix", default=".sem.sac", help="suffix of SEM sac files"
    )
    parser.add_argument(
        "--data_type", default="DISP", help="data type of SEM sac files, e.g. VEL, DISP"
    )

    return parser.parse_args()


def apply_lowpass_filter(
    trace: Trace, low_pass_freq: float, low_pass_order: int
) -> None:
    """
    Apply low-pass filter to a trace using a Butterworth filter.
    
    Args:
        trace: ObsPy trace to filter
        low_pass_freq: Corner frequency for the low-pass filter
        low_pass_order: Order of the Butterworth filter
    """
    fs = trace.stats.sampling_rate
    
    # Calculate padding and FFT length
    npad = int(2.0 * fs / low_pass_freq)
    nfft = scipy.fft.next_fast_len(trace.stats.npts + npad)
    freqs = np.fft.rfftfreq(nfft, d=1 / fs)
    
    # Apply low-pass filter in frequency domain
    lp_sos = scipy.signal.butter(
        low_pass_order, low_pass_freq, "lowpass", fs=fs, output="sos"
    )
    _, h_lp = scipy.signal.freqz_sos(lp_sos, worN=freqs, fs=fs)
    
    # Transform to frequency domain, apply filter, transform back
    sig_spectrum = np.fft.rfft(trace.data, nfft)
    sig_spectrum *= abs(h_lp)
    filtered_data = np.fft.irfft(sig_spectrum, nfft)[: trace.stats.npts]
    
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


def load_station_data(station_file: str) -> pd.DataFrame:
    """
    Load station information from CSV file.
    
    Args:
        station_file: Path to station file
        
    Returns:
        DataFrame with station information
    """
    return pd.read_csv(station_file, sep=r"\s+", header=None, comment="#")


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
        
    Returns:
        ObsPy Stream object with processed traces, or None if failed
    """
    st = Stream()
    
    # Load all channel traces
    for cha in channels:
        sac_file = os.path.join(sac_dir, f"{net}.{sta}.{loc}.{cha}{sac_suffix}")
        if not os.path.exists(sac_file):
            warnings.warn(f"[WARN] {sac_file} does not exist, skipping station {net}.{sta}.")
            return None
            
        try:
            tr = obspy.read(sac_file)
            st.extend(tr)
        except Exception as e:
            warnings.warn(f"[WARN] Failed to read {sac_file}: {e}, skipping station {net}.{sta}.")
            return None

    # Determine common time window for resampling
    max_starttime = min(tr.stats.starttime for tr in st)
    min_endtime = max(tr.stats.endtime for tr in st)
    resample_starttime = max_starttime
    resample_npts = int((min_endtime - resample_starttime) * resample_fs)

    # Process each trace: filter and resample
    for tr in st:
        apply_lowpass_filter(tr, low_pass_freq, low_pass_order)
        interpolate_trace(tr, resample_fs, resample_starttime, resample_npts)

    return st


def create_h5_dataset(
    h5f: tables.File,
    group_name: str,
    array_name: str,
    data: np.ndarray,
    attrs: dict,
    filters: tables.Filters
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
    ca = h5f.create_carray(h5f.root[group_name], array_name, atom, data.shape, filters=filters)
    ca[:] = data.astype(np.float32)
    
    # Set attributes
    for key, value in attrs.items():
        setattr(ca.attrs, key, value)
        
    return ca


def sac2h5(
    sac_dir: str,
    station_file: str,
    resample_fs: float = 2.0,
    low_pass_freq: float = 0.5,
    low_pass_order: int = 15,
    channels: List[str] = None,
    sac_suffix: str = ".sem.sac",
    data_type: str = "DISP",
    h5_file: str = "out.h5",
):
    """
    Convert SAC files to HDF5 format with filtering and resampling.
    
    Args:
        sac_dir: Directory containing SAC files
        station_file: File with station information
        resample_fs: Target sampling frequency
        low_pass_freq: Low-pass filter corner frequency
        low_pass_order: Low-pass filter order
        channels: List of channel names to process
        sac_suffix: Suffix for SAC files
        data_type: Data type (DISP, VEL, etc.)
        h5_file: Output HDF5 file path
    """
    if channels is None:
        channels = ["BXE", "BXN", "BXZ"]
    
    # Open HDF5 file
    with tables.open_file(h5_file, mode="w") as h5f:
        groot = h5f.root
        stations = load_station_data(station_file)

        for i, station in stations.iterrows():
            net, sta, stla, stlo, stel, stdp = list(station)

            print(f"[INFO] processing {net}.{sta}")

            # Extract location code from station name
            fields = sta.split(".")
            if len(fields) == 2:
                sta, loc = fields
            else:
                loc = ""

            # Process all traces for this station
            st = process_station_traces(
                sac_dir, net, sta, loc, channels, sac_suffix,
                resample_fs, low_pass_freq, low_pass_order
            )
            
            if st is None:
                continue  # Skip if traces couldn't be processed

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
                (tr.stats.channel, #.encode('utf-8'), 
                 float(tr.stats.sac["cmpaz"]), 
                 float(tr.stats.sac["cmpinc"] - 90))
                for tr in st
            ]
            
            # Create dataset and set attributes
            array_name = f"DATA_{data_type}"
            attrs = {
                
                "channels": np.array(channel_info, dtype=dtype),
                "starttime": st[0].stats.starttime,
                "sampling_rate": resample_fs,
                "npts": npts,
                "filter": [{"type": "lowpass", "N": low_pass_order, "Wn": low_pass_freq}],
                "response_corner_frequency": [],
                "type": data_type
            }
            
            create_h5_dataset(h5f, station_name, array_name, waveform_data, attrs, compression_filter)


def main():
    """Main entry point."""
    args = parse_arguments()
    print(f"Arguments: {args}")

    sac2h5(
        args.sac_dir,
        args.station,
        resample_fs=args.resample_fs,
        low_pass_freq=args.low_pass_freq,
        low_pass_order=args.low_pass_order,
        channels=args.channels,
        sac_suffix=args.sac_suffix,
        data_type=args.data_type,
        h5_file=args.out_h5,
    )


if __name__ == "__main__":
    main()