"""
Test if events are inside the SEM mesh domain.

This script reads a SPECFEM3D_GLOBE parameter file and an events CSV file,
converts event coordinates to mesh xi-eta coordinates, and outputs the results.
"""

import argparse
import pandas as pd

from meshfem3d_utils import sem_latlon2xieta


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
        "events_csvfile",
        type=str,
        help="Events CSV file with pipe-separated columns: gcmtID|time|lat|lon|depth|Mw|description",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="events_with_xi_eta.csv",
        help="Output CSV file path (default: events_with_xi_eta.csv), pipe-separated",
    )
    return parser.parse_args()


def parse_sem_parameters(parfile: str) -> dict:
    """Parse SPECFEM3D_GLOBE Par_file to a dictionary."""
    sem_params = pd.read_csv(
        parfile,
        delimiter=r"\s*=\s*",
        header=None,
        comment="#",
        names=["key", "value"],
        dtype=dict(key=object, value=object),
        index_col=["key"],
        engine="python",
    ).to_dict()["value"]

    return sem_params
    

def compute_xi_eta(events_df: pd.DataFrame, sem_params: dict) -> pd.DataFrame:
    """Compute xi-eta coordinates for all events."""
    # Convert scientific notation (e.g., "1.0d2" -> "1.0e2") and extract mesh parameters
    mesh_central_lat = float(sem_params['CENTER_LATITUDE_IN_DEGREES'].lower().replace('d', 'e'))
    mesh_central_lon = float(sem_params['CENTER_LONGITUDE_IN_DEGREES'].lower().replace('d', 'e'))
    mesh_width_xi = float(sem_params['ANGULAR_WIDTH_XI_IN_DEGREES'].lower().replace('d', 'e'))
    mesh_width_eta = float(sem_params['ANGULAR_WIDTH_ETA_IN_DEGREES'].lower().replace('d', 'e'))
    mesh_gamma_rot = float(sem_params['GAMMA_ROTATION_AZIMUTH'].lower().replace('d', 'e'))
    
    def _get_xi_eta(row):
        xi, eta = sem_latlon2xieta(
            mesh_central_lat,
            mesh_central_lon,
            mesh_gamma_rot,
            row["lat"],
            row["lon"],
        )
        dist_from_edge = min(0.5*mesh_width_xi - abs(xi), 0.5*mesh_width_eta - abs(eta))
        return pd.Series([xi, eta, dist_from_edge], index=["xi", "eta", "dist_from_edge"])
    
    return events_df.apply(_get_xi_eta, axis=1)


if __name__ == "__main__":
    args = parse_arguments()
    
    # Parse mesh parameters
    sem_params = parse_sem_parameters(args.sem_parfile)
    
    # Read events
    events_df = pd.read_csv(args.events_csvfile, sep="|", header=None)
    events_df.columns = ["gcmtID", "time", "lat", "lon", "depth", "mag", "description"]
    
    # Compute xi-eta coordinates
    df_xi_eta = compute_xi_eta(events_df, sem_params)
    
    # Save results
    result_df = pd.concat([events_df, df_xi_eta], axis=1)
    result_df.to_csv(args.output, index=False, sep="|")
    print(f"Results saved to {args.output}")
