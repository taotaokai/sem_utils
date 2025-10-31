#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert CMTSOLUTION file from geographic coordinates to ECEF coordinates.

This script converts a CMTSOLUTION file from geographic coordinates (latitude, longitude, depth)
to Earth-Centered Earth-Fixed (ECEF) coordinates. It also transforms the moment tensor from
the radial-transverse-phi (r,theta,phi) coordinate system to the Cartesian (x,y,z) system.

The output includes:
- Event location in ECEF coordinates (meters)
- Moment tensor in ECEF Cartesian coordinates
- Converted time parameters (t0=0, tau parameter)
"""

import os
import argparse
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from netCDF4 import Dataset
from obspy import UTCDateTime
from pyproj import Transformer
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_cmt_file(cmt_file):
    """
    Parse the CMTSOLUTION file and extract parameters.

    Args:
        cmt_file (str): Path to the CMTSOLUTION file

    Returns:
        dict: Dictionary containing parsed parameters
    """
    with open(cmt_file, "r") as f:
        lines = [x for x in f.readlines() if not x.startswith("#")]

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
    tau = float(lines[3][1]) / 1.628  # mimic triangle with gaussian
    lat = float(lines[4][1])
    lon = float(lines[5][1])
    dep = float(lines[6][1]) * 1000.0  # Convert from km to meters

    # Moment tensor components (convert from dyn*cm to N*m by multiplying with 1e-7)
    m11 = float(lines[7][1])
    m22 = float(lines[8][1])
    m33 = float(lines[9][1])
    m12 = float(lines[10][1])
    m13 = float(lines[11][1])
    m23 = float(lines[12][1])

    mt_rtp = np.array([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]]) * 1.0e-7

    return {
        "header": header,
        "year": year,
        "month": month,
        "day": day,
        "hour": hour,
        "minute": minute,
        "second": second,
        "event_id": event_id,
        "time_shift": time_shift,
        "tau": tau,
        "lat": lat,
        "lon": lon,
        "dep": dep,
        "mt_rtp": mt_rtp,
    }


def load_topography(topo_ncfile):
    """
    Load topography data from NetCDF file.

    Args:
        topo_ncfile (str): Path to the topography NetCDF file

    Returns:
        RegularGridInterpolator: Interpolator for topography data or None if file not provided
    """
    if not topo_ncfile:
        return None

    if not os.path.exists(topo_ncfile):
        logger.error(f"Topography file {topo_ncfile} does not exist")
        return None

    try:
        topo_ds = Dataset(topo_ncfile, "r")
        topo_heights = np.array(topo_ds.variables["z"])
        topo_lons = np.array(topo_ds.variables["lon"])
        topo_lats = np.array(topo_ds.variables["lat"])
        interpolator = RegularGridInterpolator((topo_lats, topo_lons), topo_heights)
        topo_ds.close()
        return interpolator
    except Exception as e:
        logger.error(f"Error loading topography file {topo_ncfile}: {e}")
        return None


def lla_to_ecef(lat, lon, height):
    """
    Convert latitude, longitude, altitude to ECEF coordinates.

    Args:
        lat (float): Latitude in degrees
        lon (float): Longitude in degrees
        height (float): Height above ellipsoid in meters

    Returns:
        tuple: ECEF coordinates (x, y, z) in meters
    """
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:4978")  # WGS84 to ECEF
    return transformer.transform(lat, lon, height)


def calculate_ecef_coordinates(lat, lon, dep, topo_interp=None):
    """
    Calculate ECEF coordinates accounting for topography.

    Args:
        lat (float): Latitude in degrees
        lon (float): Longitude in degrees
        dep (float): Depth in meters (positive down)
        topo_interp (RegularGridInterpolator): Topography interpolator or None

    Returns:
        tuple: (x, y, z) ECEF coordinates in meters
    """
    # Get WGS84 height
    topo = 0.0
    if topo_interp:
        try:
            topo = topo_interp((lat, lon))  # Get local topography
        except Exception as e:
            logger.warning(f"Failed to interpolate topography at {lat}, {lon}: {e}")

    height = topo - dep  # WGS84 ellipsoidal height of the earthquake
    logger.info(f"Coordinates: lon={lon}, lat={lat}, topo={topo}, height={height}")

    # Convert from LLA to ECEF (meters)
    return lla_to_ecef(lat, lon, height)


def compute_coordinate_transformation(x, y, z):
    """
    Compute coordinate transformation matrix from (r,theta,phi) to (x,y,z).

    Args:
        x, y, z (float): ECEF coordinates in meters

    Returns:
        np.array: 3x3 transformation matrix
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    sthe = np.sin(theta)
    cthe = np.cos(theta)
    sphi = np.sin(phi)
    cphi = np.cos(phi)

    return np.array(
        [
            [sthe * cphi, cthe * cphi, -sphi],
            [sthe * sphi, cthe * sphi, cphi],
            [cthe, -sthe, 0.0],
        ]
    )


def transform_moment_tensor(mt_rtp, transformation_matrix):
    """
    Transform moment tensor from radial-transverse-phi to Cartesian coordinates.

    Args:
        mt_rtp (np.array): 3x3 moment tensor in radial-transverse-phi coordinates
        transformation_matrix (np.array): 3x3 coordinate transformation matrix

    Returns:
        np.array: 3x3 moment tensor in Cartesian coordinates
    """
    return np.dot(np.dot(transformation_matrix, mt_rtp), transformation_matrix.T)


def calculate_centroid_time(params):
    """
    Calculate the centroid time by applying time shift to origin time.

    Args:
        params (dict): Dictionary containing time parameters

    Returns:
        UTCDateTime: Centroid time
    """
    isotime = f"{params['year']}-{params['month']}-{params['day']}T{params['hour']}:{params['minute']}:{params['second']}Z"
    return UTCDateTime(isotime) + params["time_shift"]


def format_header_line(params, centroid_time):
    """
    Format the header line with the centroid time.

    Args:
        params (dict): Dictionary containing original header parameters
        centroid_time (UTCDateTime): Calculated centroid time

    Returns:
        str: Formatted header line
    """
    header = params["header"]
    header[1] = f"{centroid_time.year:04d}"
    header[2] = f"{centroid_time.month:02d}"
    header[3] = f"{centroid_time.day:02d}"
    header[4] = f"{centroid_time.hour:02d}"
    header[5] = f"{centroid_time.minute:02d}"
    header[6] = (
        f"{centroid_time.second:02d}.{round(centroid_time.microsecond * 1e-3):03d}"
    )
    return " ".join(header)


def write_output_file(out_file, header_line, params, x, y, z, mt_xyz, tau):
    """
    Write the converted CMTSOLUTION file in ECEF coordinates.

    Args:
        out_file (str): Output file path
        header_line (str): Formatted header line
        params (dict): Dictionary containing event parameters
        x, y, z (float): ECEF coordinates in meters
        mt_xyz (np.array): Moment tensor in ECEF coordinates
        tau (float): Tau parameter
    """
    try:
        with open(out_file, "w") as fp:
            fp.write(f"{header_line}\n")
            fp.write(f"{'event_name:':<18} {params['event_id']}\n")
            fp.write(f"{'t0(s):':<18} {0.0:+15.8E}\n")
            fp.write(f"{'tau(s):':<18} {tau:+15.8E}\n")
            fp.write(f"{'x(m):':<18} {x:+15.8E}\n")
            fp.write(f"{'y(m):':<18} {y:+15.8E}\n")
            fp.write(f"{'z(m):':<18} {z:+15.8E}\n")
            fp.write(f"{'Mxx(N*m):':<18} {mt_xyz[0, 0]:+15.8E}\n")
            fp.write(f"{'Myy(N*m):':<18} {mt_xyz[1, 1]:+15.8E}\n")
            fp.write(f"{'Mzz(N*m):':<18} {mt_xyz[2, 2]:+15.8E}\n")
            fp.write(f"{'Mxy(N*m):':<18} {mt_xyz[0, 1]:+15.8E}\n")
            fp.write(f"{'Mxz(N*m):':<18} {mt_xyz[0, 2]:+15.8E}\n")
            fp.write(f"{'Myz(N*m):':<18} {mt_xyz[1, 2]:+15.8E}\n")
        logger.info(f"Successfully wrote output to {out_file}")
    except Exception as e:
        logger.error(f"Error writing output file {out_file}: {e}")
        raise


def main():
    """Main function to convert CMTSOLUTION from geographic to ECEF coordinates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cmt_file", help="Input CMTSOLUTION file")
    parser.add_argument("out_file", help="Output CMTSOLUTION file in ECEF")
    parser.add_argument(
        "--topo_ncfile",
        default=None,
        help="Topo netcdf file with lat,lon,z[lat,lon] variables. Vertical reference to WGS84 ellipsoid.",
    )

    args = parser.parse_args()
    logger.info(f"Arguments: {args}")

    # Validate input file
    if not os.path.exists(args.cmt_file):
        logger.error(f"Input CMTSOLUTION file not found: {args.cmt_file}")
        return 1

    # Load topography data if provided
    topo_interp = load_topography(args.topo_ncfile)
    if args.topo_ncfile and not topo_interp:
        logger.warning("Failed to load topography data, continuing without it")

    # Parse CMTSOLUTION file
    try:
        params = parse_cmt_file(args.cmt_file)
    except Exception as e:
        logger.error(f"Error parsing CMTSOLUTION file {args.cmt_file}: {e}")
        return 1

    # Calculate ECEF coordinates
    try:
        x, y, z = calculate_ecef_coordinates(
            params["lat"], params["lon"], params["dep"], topo_interp
        )
    except Exception as e:
        logger.error(f"Error calculating ECEF coordinates: {e}")
        return 1

    # Compute coordinate transformation matrix
    try:
        transformation_matrix = compute_coordinate_transformation(x, y, z)
    except Exception as e:
        logger.error(f"Error computing coordinate transformation: {e}")
        return 1

    # Transform moment tensor
    try:
        mt_xyz = transform_moment_tensor(params["mt_rtp"], transformation_matrix)
    except Exception as e:
        logger.error(f"Error transforming moment tensor: {e}")
        return 1

    # Calculate centroid time
    try:
        centroid_time = calculate_centroid_time(params)
    except Exception as e:
        logger.error(f"Error calculating centroid time: {e}")
        return 1

    # Format header with centroid time
    header_line = format_header_line(params, centroid_time)

    # Write output file
    try:
        write_output_file(
            args.out_file, header_line, params, x, y, z, mt_xyz, params["tau"]
        )
    except Exception as e:
        logger.error(f"Error writing output file: {e}")
        return 1

    logger.info("Conversion completed successfully")
    return 0


if __name__ == "__main__":
    exit(main())
