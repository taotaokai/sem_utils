import argparse
import numpy as np
import pandas as pd
import pyproj
from obspy.taup import TauPyModel

"""Parse command line arguments."""
parser = argparse.ArgumentParser(
    description="decluster events by distances and time interval"
)

parser.add_argument("event_list", help="events list to filter")
parser.add_argument("gcmt_list", help="global gcmt list")
parser.add_argument("out_list", help="selected list")
parser.add_argument("--central_lat", type=float, required=True, help="central latitude")
parser.add_argument(
    "--central_lon", type=float, required=True, help="central longitude"
)
parser.add_argument(
    "--time_before",
    type=float,
    default=6000,
    help="time before event when no other event occurs",
)
parser.add_argument(
    "--time_after",
    type=float,
    default=1800,
    help="time after event when no other event occurs",
)
parser.add_argument(
    "--amplitude_ratio",
    type=float,
    default=0.05,
    help="minimum amplitue ratio threshold to ignore effect of overlapping events",
)
args = parser.parse_args()

gps2ecef = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")  # GPS to ECEF
# ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")  # ECEF to GPS
geod = pyproj.Geod(ellps="WGS84")
R_EARTH = (geod.a + geod.b) / 2

taup = TauPyModel(model="ak135")

df_gcmt = pd.read_csv(
    args.gcmt_list,
    sep="|",
    header=None,
    names=["gcmtID", "time", "lat", "lon", "depth", "mag", "description"],
    comment="#",
)
df_gcmt["time"] = pd.to_datetime(df_gcmt["time"])

df = pd.read_csv(
    args.event_list,
    sep="|",
    header=None,
    names=["gcmtID", "time", "lat", "lon", "depth", "mag", "description"],
)
df["time"] = pd.to_datetime(df["time"])
df.sort_values("mag", inplace=True, ascending=False)

df_sel = df.iloc[0:0].copy()  # empty dataframe with same columns and dtypes

for i, row in df.iterrows():
    print("============ process event\n", row)
    event_time = row["time"]
    event_mag = row["mag"]
    event_lat = row["lat"]
    event_lon = row["lon"]
    event_dep = row["depth"]

    _, _, dist = geod.inv(event_lon, event_lat, args.central_lon, args.central_lat)
    dist = np.rad2deg(dist / R_EARTH)

    # arrs = taup.get_travel_times(event_lat, dist, "ttp")
    # min_ttp = np.min([arr.time for arr in arrs])
    # min_ttsurf = 30 * dist

    t0 = event_time - pd.Timedelta(args.time_before, unit="seconds")
    t1 = event_time + pd.Timedelta(args.time_after, unit="seconds")

    ind = (df_gcmt["time"] >= t0) & (df_gcmt["time"] <= t1) & (df_gcmt["gcmtID"] != row["gcmtID"])

    df_overlap = df_gcmt[ind]

    if len(df_overlap) == 0:
        df_sel.loc[len(df_sel), :] = row[:]
        continue

    print("------ overlapping events\n", df_overlap)

    n = len(df_overlap)
    _, _, other_dists = geod.inv(
        df_overlap["lon"].values, df_overlap["lat"].values, [args.central_lon]*n, [args.central_lat]*n
    )
    other_dists = np.rad2deg(np.array(other_dists) / R_EARTH)
    other_mags = df_overlap["mag"].values

    amplitude_ratio = 10 ** (other_mags - event_mag) * (other_dists / dist) ** (-1.66)

    print(f"{dist=}")
    print(f"{other_dists=}")
    print(f"{amplitude_ratio=}")

    if np.any(amplitude_ratio > args.amplitude_ratio):
        print("!!!!!!! EVENT IGNORED !!!!!!!")
    else:
        df_sel.loc[len(df_sel), :] = row[:]

df_sel.to_csv(args.out_list, sep="|", header=None, index=False)