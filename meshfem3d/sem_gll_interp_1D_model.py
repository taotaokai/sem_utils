#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Interpolate 1D velocity model to SEM GLL points.

This script interpolates a 1D velocity model (defined by depth profiles)
onto the Gauss-Lobatto-Legendre (GLL) points of a spectral element mesh.
The interpolation is performed for each processor's mesh partition.

Author: [Author Name]
Date: [Date]
"""

import argparse
import sys
import os
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np
import pyproj

from meshfem3d_constants import R_EARTH
from meshfem3d_utils import sem_mesh_read, write_gll_file


def parse_command_line_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for the 1D model interpolation script.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Interpolate 1D velocity model to SEM GLL points",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--model_csv",
        required=True,
        help="1-D model CSV file with depth_km, v1, v2 columns",
    )
    parser.add_argument(
        "--nproc", required=True, type=int, help="Number of processors/processes"
    )
    parser.add_argument(
        "--mesh_dir", default="DATABASES_MPI", help="Directory containing mesh files"
    )
    parser.add_argument(
        "--model_names",
        nargs="+",
        required=True,
        help="Parameter names in the model (e.g., vp vs rho)",
    )
    parser.add_argument(
        "--out_dir",
        default="interp_model",
        help="Output directory for interpolated models",
    )

    return parser.parse_args()


def read_model_file(
    model_csv: str, model_names: List[str]
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Read the 1D model CSV file and extract depth and parameter values.

    Args:
        model_csv (str): Path to the model CSV file
        model_names (List[str]): List of parameter names to extract

    Returns:
        Tuple[np.ndarray, Dict[str, np.ndarray]]: Depth array and dictionary of parameter values

    Raises:
        FileNotFoundError: If model CSV file doesn't exist
        ValueError: If required columns are missing
    """
    if not os.path.exists(model_csv):
        raise FileNotFoundError(f"Model file not found: {model_csv}")

    try:
        model_df = pd.read_csv(model_csv)
    except Exception as e:
        raise ValueError(f"Error reading model file: {e}")

    # Check if required columns exist
    required_columns = ["depth_km"] + model_names
    missing_columns = [col for col in required_columns if col not in model_df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in model file: {missing_columns}")

    # Extract depth and parameter values
    model_depths = np.array(model_df["depth_km"])
    model_values = {}

    for param_name in model_names:
        model_values[param_name] = np.array(model_df[param_name])

    return model_depths, model_values


def process_mesh_slice(
    iproc: int,
    mesh_dir: str,
    out_dir: str,
    model_depths: np.ndarray,
    model_values: Dict[str, np.ndarray],
    model_names: List[str],
    ecef2gps_transformer,
) -> None:
    """
    Process mesh data for a single processor and interpolate model values.

    Args:
        iproc (int): Processor number
        mesh_dir (str): Directory containing mesh files
        out_dir (str): Output directory
        model_depths (np.ndarray): Array of model depths
        model_values (Dict[str, np.ndarray]): Dictionary of model parameter values
        model_names (List[str]): List of parameter names
        ecef2gps_transformer: PyProj transformer for coordinate conversion
    """
    print(f"Processing slice {iproc:06d}")

    # Read target SEM mesh
    mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")

    try:
        mesh = sem_mesh_read(mesh_file)
    except Exception as e:
        print(f"Warning: Could not read mesh file {mesh_file}: {e}")
        return

    ibool = mesh["ibool"]
    xyz_glob = mesh["xyz_glob"] * R_EARTH

    # Convert ECEF coordinates to geographic coordinates
    lat, lon, alt = ecef2gps_transformer.transform(
        xyz_glob[:, 0], xyz_glob[:, 1], xyz_glob[:, 2]
    )

    # Calculate depth from altitude (convert meter to km)
    depth_glob = -alt / 1000.0

    # Interpolate model values for each parameter
    for param_name in model_names:
        try:
            model_interp = np.interp(depth_glob, model_depths, model_values[param_name])

            # Save interpolated model to GLL file
            write_gll_file(
                out_dir,
                param_name,
                iproc,
                model_interp[ibool],
                overwrite=True,
            )

        except Exception as e:
            print(
                f"Warning: Error processing parameter {param_name} for processor {iproc}: {e}"
            )


def main():
    """
    Main function to orchestrate the 1D model interpolation process.
    """
    try:
        # Parse command line arguments
        args = parse_command_line_arguments()
        print("Running with arguments:", args)

        # Create output directory if it doesn't exist
        os.makedirs(args.out_dir, exist_ok=True)

        # Initialize coordinate transformer (ECEF to geographic)
        ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")

        # Read model file
        print("Reading model file...")
        model_depths, model_values = read_model_file(args.model_csv, args.model_names)
        print(f"Model loaded with {len(model_depths)} depth levels")

        # Process each processor's mesh
        print("Starting interpolation process...")
        for iproc in range(args.nproc):
            process_mesh_slice(
                iproc,
                args.mesh_dir,
                args.out_dir,
                model_depths,
                model_values,
                args.model_names,
                ecef2gps,
            )

        print("Interpolation completed successfully!")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
