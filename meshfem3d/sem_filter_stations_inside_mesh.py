"""
Test if events are inside the SEM mesh domain.

This script reads a SPECFEM3D_GLOBE parameter file and an events CSV file,
converts event coordinates to mesh xi-eta coordinates, and outputs the results.
"""

import argparse

import numpy as np
import pandas as pd
from meshfem3d_utils import read_sem_parfile, sem_latlon2xieta


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert event coordinates to SEM mesh xi-eta coordinates and test if they are inside the mesh domain."
    )
    parser.add_argument(
        "sem_parfile",
        type=str,
        help="SPECFEM3D_GLOBE parameter file (e.g., Par_file)",
    )
    parser.add_argument(
        "station_file",
        type=str,
        help="station file with lat,lon",
    )
    parser.add_argument(
        "--min_dist",
        required=True,
        type=float,
        default=1.0,
        help="minimum distance from mesh boundaries (in degree)",
    )
    parser.add_argument(
        "--lat_pos",
        required=True,
        type=int,
        help="the column position of latitudes",
    )
    parser.add_argument(
        "--lon_pos",
        required=True,
        type=int,
        help="the column position of latitudes",
    )
    parser.add_argument(
        "--sep",
        type=str,
        default=r"\s+",
        help="column seperator for the input/output files",
    )
    parser.add_argument(
        "--header",
        type=int,
        default=None,
        help="number of header lines",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="filtered.csv",
        help="Output CSV file path",
    )
    parser.add_argument(
        "--out_delimiter",
        type=str,
        default=" ",
        help="delimiter of the output file",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Parse mesh parameters
    sem_params = read_sem_parfile(args.sem_parfile)
    mesh_central_lat = sem_params["CENTER_LATITUDE_IN_DEGREES"]
    mesh_central_lon = sem_params["CENTER_LONGITUDE_IN_DEGREES"]
    mesh_width_xi = sem_params["ANGULAR_WIDTH_XI_IN_DEGREES"]
    mesh_width_eta = sem_params["ANGULAR_WIDTH_ETA_IN_DEGREES"]
    mesh_gamma_rot = sem_params["GAMMA_ROTATION_AZIMUTH"]

    # Read input station file
    df = pd.read_csv(args.station_file, sep=args.sep, header=args.header)

    # Compute xi-eta coordinates
    lats = df.iloc[:, args.lat_pos]
    lons = df.iloc[:, args.lon_pos]
    xi, eta = sem_latlon2xieta(
        mesh_central_lat, mesh_central_lon, mesh_gamma_rot, lats, lons
    )
    dist_from_edge = np.minimum(0.5 * mesh_width_xi - np.abs(xi), 0.5 * mesh_width_eta - np.abs(eta))

    df_filter = df.iloc[dist_from_edge >= args.min_dist]

    # Save results
    df_filter.to_csv(args.output, index=False, header=False, sep=args.out_delimiter)
    print(f"Results saved to {args.output}")
