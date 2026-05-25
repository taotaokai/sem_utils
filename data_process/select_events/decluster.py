import argparse
import pandas as pd
import pyproj

"""Parse command line arguments."""
parser = argparse.ArgumentParser(
    description="decluster events by distances and time interval"
)

parser.add_argument("event_list", help="events list")
parser.add_argument("out_list", help="output list of selected events")
parser.add_argument("--min_dist_in_km", type=int, default=100, help="minimum distance in km")
parser.add_argument("--min_time_diff_in_days", type=int, default=180, help="minimum time interval in days")
args = parser.parse_args()

gps2ecef = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:4978")  # GPS to ECEF
# ecef2gps = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326")  # ECEF to GPS

# event_list = "events_sel_2.lst"
# min_dist_in_km = 100
# min_time_diff_in_sec = 6 * 30 * 24 * 3600 # 6 month 

df = pd.read_csv(
    args.event_list,
    sep="|",
    header=None,
    names=["gmtID", "time", "lat", "lon", "depth", "mag", "description"],
)
df["time"] = pd.to_datetime(df["time"])
df.sort_values("mag", inplace=True, ascending=False)

lats = df["lat"].values
lons = df["lon"].values
alts = -1000.0 * df["depth"].values

xx, yy, zz = gps2ecef.transform(lats, lons, alts)

df["x"] = xx
df["y"] = yy
df["z"] = zz

df_sel = df.iloc[0:0].copy() # empty dataframe with same columns and dtypes

for i, row in df.iterrows():
    print("============ process event\n", row)
    event_time = row["time"]
    event_mag = row["mag"]
    event_lat = row["lat"]
    event_lon = row["lon"]

    if event_mag > 7.0:
        continue    

    if len(df_sel) == 0:
        df_sel.loc[len(df_sel), :] = row[:] # NOTE must use [,:] = row[:] to keep dtypes
        continue

    # distance
    dx = df_sel["x"] - row["x"]
    dy = df_sel["y"] - row["y"]
    dz = df_sel["z"] - row["z"]

    dist = (dx**2 + dy**2 + dz**2)**0.5 / 1000.0 # km

    df1 = df_sel[dist <= args.min_dist_in_km]
    print("----- nearby events: \n", df1)

    if len(df1) == 0:
        df_sel.loc[len(df_sel)] = row
        continue

    dt = df1["time"] - event_time
    ind = abs(dt) < pd.Timedelta(args.min_time_diff_in_days, unit="day")
    if any(ind):
        print(f"----- dumplicated events:\n{df1[ind]}")
        print("!!!!!!!! DROP EVENT")
        continue

    df_sel.loc[len(df_sel)] = row

df_sel.drop(columns=["x", "y", "z"], inplace=True)
df_sel['time'] = df_sel['time'].dt.strftime('%Y-%m-%dT%H:%M:%S.%f')

df_sel.to_csv(args.out_list, sep="|", header=None, index=False)