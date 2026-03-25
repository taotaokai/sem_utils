import os
import numpy as np
import tables
import scipy
import yaml
import argparse
import logging

from obspy import read_inventory, Inventory, read_events, geodetics
from obspy.clients.filesystem.sds import Client
from obspy.taup import TauPyModel


""" data structure

/root
    |- attrs: event_id="C202006211907A", origin_time='2020-06-21T19:07:00.000000', ...
    |
    |- BE_MEM
    |   |- attrs: network='BE', station='MEM', location='00',
    |   |         stla,stlo,stel,stdp
    |   |
    |   |- DATA_RAW[0:3, 0:npts] (carray)
    |       |- attrs: channels=[('HHE',90,0), ('HHN',0,0), ('HHZ',0,-90)],
    |       |         filter=[{'type': 'lowpass', 'N': 15, 'Wn': 0.5}, 
    |       |                 {'type': 'highpass', 'N': 2, 'Wn': 0.033}]
    |       |         sampling_rate=10.0,
    |       |         npts=18000,
    |       |         unit='count',
    |       |         type='raw'
    |       |         starttime='2020-03-21T00:49:51.008393',
    |
    |- (other stations) 

NOTE: time precision is microsecond

"""

args = argparse.ArgumentParser()
args.add_argument("config_yaml", help="configuration yaml file")
args.add_argument("CMT_file", help="CMT file")
args.add_argument("out_h5_file", help="output h5 file")
args.add_argument("out_channel_file", help="output channel file")
args.add_argument("--log_file", default=None, help="log file")
args.add_argument("--quiet", action="store_true", help="suppress output to console")
args = args.parse_args()

config_yaml = args.config_yaml
CMT_file = args.CMT_file
out_h5_file = args.out_h5_file
out_channel_file = args.out_channel_file
log_file = args.log_file

# Create a custom logger
logger = logging.getLogger("SDS_cut_event_h5")  # Use __name__ for best practices
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d_%H:%M:%S"
)
# Create handlers
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
if args.quiet:
    console_handler.setLevel(logging.ERROR)
else:
    console_handler.setLevel(logging.INFO)
# Add handlers to the logger
logger.addHandler(console_handler)

if log_file:
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

# print("==================================\n")
# print(f"START: {datetime.datetime.now()}\n")
logger.info("START >>>>>>>>>>>>>>>>>>")
logger.info(f"{args}")

if os.path.dirname(out_h5_file):
    os.makedirs(os.path.dirname(out_h5_file), exist_ok=True)
if os.path.dirname(out_channel_file):
    os.makedirs(os.path.dirname(out_channel_file), exist_ok=True)

# h5 file handle
h5f = tables.open_file(out_h5_file, mode="w")

# event
event_data = read_events(CMT_file)[0]
event_id = event_data.event_descriptions[0].text
event_origin = event_data.preferred_origin()

logger.info(f"{event_origin}")

# config
with open(config_yaml, "r") as file:
    config = yaml.safe_load(file)

config_zcomp_pattern = config["zcomp_pattern"]

config_sampling_rate = config["resample"]["sampling_rate"]
config_lowpass_N = config["resample"]["filter"]["N"]
config_lowpass_Wn = config["resample"]["filter"]["Wn"]

time_before_first_arrival = config["time_window"]["before_first_arrival"]
if not isinstance(time_before_first_arrival, list):
    time_before_first_arrival = [
        time_before_first_arrival,
    ]
# time_after_first_arrival = config['data']['time_window']['after_first_arrival_seconds']
time_after_origin = config["time_window"]["after_origin_time"]
if not isinstance(time_after_origin, list):
    time_after_origin = [
        time_after_origin,
    ]

# for merging traces
# misalignment_threshold = config['misalignment_threshold']

stations_path = config["stations_path"]  # station xml
SDS_root = config["SDS_root"]  # mseed storage

groot = h5f.root

# save configuration
groot._v_attrs["config"] = config

# set attributes for event
groot._v_attrs["event"] = event_data
groot._v_attrs["event_id"] = event_id
groot._v_attrs["origin_time"] = event_origin.time
groot._v_attrs["location"] = {
    "longitude": event_origin.longitude,
    "latitude": event_origin.latitude,
    "depth_km": event_origin.depth / 1000.0,
}

# taup model
taup_model = TauPyModel(model="ak135")

# SDS client
client = Client(SDS_root)

t0 = event_origin.time - max(time_before_first_arrival)
t1 = event_origin.time + max(time_after_origin)

# query all Z components
st = client.get_waveforms("*", "*", "*", config_zcomp_pattern, t0, t1)

stations = [(tr.stats.network, tr.stats.station) for tr in st]
stations = sorted(set(stations))  # remove duplicates

nslc_list = [
    (tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel)
    for tr in st
]
nslc_list = sorted(set(nslc_list))  # remove duplicates

# in case multiple location/band, keep the one with largest coverage
nslc_filtered = []
for net, sta in stations:
    sxml_file = os.path.join(stations_path, f"{net}.{sta}.xml")
    try:
        inv = read_inventory(sxml_file)
    except Exception as err:
        logger.warning(f"failed to read {sxml_file}, skip. (Error: {err})")
        continue
    # inv = inv.select(network=net, station=sta, time=event_origin.time)
    # if not inv:
    #     print(f"[WARN] fail to select inventory for {net}.{sta} on {event_origin.time}")
    st1 = st.select(net, sta)
    nslc_list = [
        (tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel)
        for tr in st1
    ]
    nslc_list = sorted(set(nslc_list))  # remove duplicates
    nslc_use = []
    for nslc in nslc_list:
        net, sta, loc, cha = nslc
        inv1 = inv.select(
            network=net, station=sta, location=loc, channel=cha, time=event_origin.time
        )
        if not inv1:
            logger.warning(
                f"failed to select inventory for {nslc} on {event_origin.time}"
            )
            # print(f"[WARN] fail to select inventory for {nslc} on {event_origin.time}")
            continue
        duration = sum(
            [tr.stats.endtime - tr.stats.starttime for tr in st1.select(*nslc)]
        )
        coverage = duration / (t1 - t0)
        nslc_use.append((nslc, coverage))
    if nslc_use:
        nslc_use = sorted(nslc_use, key=lambda x: x[1])  # sort by data coverage
        nslc_filtered.append(
            nslc_use[-1][0]
        )  # only keep the trace with largest coverage for each station


# process
inv_used = Inventory()
for network, station, location, channel in nslc_filtered:
    station_id = ".".join([network, station, location, channel[:-1] + "?"])
    # print(f"[INFO] {station_id}\n")
    logger.info(f"====== {station_id}")

    # read StationXML file
    sxml_file = os.path.join(stations_path, f"{network}.{station}.xml")
    try:
        inv = read_inventory(sxml_file)
    except Exception as err:
        logger.warning(f"fail to read {sxml_file}, station ignored. ({station_id})")
        logger.debug(f"Error message: {err}")
        # print(
        #     f"{err}\n[WARN] fail to read {sxml_file}, station ignored. ({station_id})\n"
        # )
        continue

    inv = inv.select(
        network=network,
        station=station,
        location=location,
        channel=channel[:-1] + "?",
        time=event_origin.time,
    )
    if not inv:
        logger.warning(
            f"no station metadata on {event_origin.time}), station ignored. ({station_id})"
        )
        # print(
        #     f"[WARN] no station metadata on {event_origin.time}), station ignored. ({station_id})\n"
        # )
        continue

    # get channel coordinates
    try:
        seed_id = ".".join([network, station, location, channel])
        station_coords = inv.get_coordinates(seed_id, event_origin.time)
    except Exception as err:
        logger.warning(
            f"cannot find station metadata for {seed_id}, station ignored. ({station_id})"
        )
        logger.debug(f"Error message: {err}")
        # print(
        #     f"{err}\n[WARN] cannot find station metadata for {seed_id}, station ignored. ({station_id})\n"
        # )
        continue

    # calculate event-station distance
    dist, az, baz = geodetics.gps2dist_azimuth(
        event_origin.latitude,
        event_origin.longitude,
        station_coords["latitude"],
        station_coords["longitude"],
    )
    dist_degree = geodetics.kilometer2degrees(dist / 1000.0)

    # calculate first arrival time using ak135 model
    arrivals = taup_model.get_travel_times(
        source_depth_in_km=event_origin["depth"] / 1000.0,
        distance_in_degree=dist_degree,
        phase_list=["ttp"],
    )
    first_arrival_time = event_origin.time + min([arr.time for arr in arrivals])

    # get waveforms
    t0 = first_arrival_time - max(time_before_first_arrival)
    t1 = time_window_endtime = event_origin.time + max(time_after_origin)
    try:
        pad_time = 5  # seconds
        st = client.get_waveforms(
            network,
            station,
            location,
            channel[0:2] + "?",
            t0,
            t1,
            merge=-1,
        )
    except Exception as err:
        logger.warning(f"failed to get waveforms, station ingored. ({station_id})")
        logger.warning(f"encounter an error: {err}")
        # print(
        #     f"{err}\n[WARN] failed to get waveforms, station ignored. ({station_id})\n"
        # )
        continue

    # merge adjacent traces with misalignment less than sub-sampling interval
    # NOTE merge is done in client.get_waveforms(), which defaults merge=-1
    #      default misalignment_threshold=0.01 should be fine
    # st.merge(-1, misalignment_threshold=misalignment_threshold)

    # ignore traces with gaps or overlaps
    orientation_codes = set(
        [tr.stats.channel[-1] for tr in st]
    )  # get unique orientation codes
    for orientation in orientation_codes:
        st1 = st.select(component=orientation)
        if len(st1) != 1:
            # print(f"{st1.__str__(extended=True)}")
            # print(f"[WARN] gaps/overlaps found, {st1[0].id} ignored. ({station_id})\n")
            logger.warning(f"gaps/overlaps found, {st1[0].id} ignored. ({station_id})")
            logger.debug(f"{st1.__str__(extended=True)}")
            for tr in st1:
                st.remove(tr)
    if len(st) == 0:
        logger.warning(f"no trace left, station ignored. ({station_id})")
        # print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # check time coverage
    traces_to_remove = []
    # [t0, t1]: minimum time range required
    t0 = first_arrival_time - min(time_before_first_arrival)
    t1 = event_origin.time + min(time_after_origin)
    for tr in st:
        dt = tr.stats.delta
        if tr.stats.starttime > t0:
            logger.warning(
                f"trace starttime {tr.stats.starttime} is later than required {t0}, {tr.id} ignored. ({station_id})"
            )
            logger.debug(f"{tr}")
            # print(
            #     f"{tr}\n[WARN] trace starttime {tr.stats.starttime} is later than required {t0}, {tr.id} ignored. ({station_id})\n"
            # )
            traces_to_remove.append(tr)
        elif tr.stats.endtime < t1:
            logger.warning(
                f"trace endtime {tr.stats.endtime} is earlier than required {t1}, {tr.id} ignored. ({station_id})"
            )
            logger.debug(f"{tr}")
            # print(
            #     f"{tr}\n[WARN] trace endtime {tr.stats.endtime} is earlier than required {t1}, {tr.id} ignored. ({station_id})\n"
            # )
            traces_to_remove.append(tr)
    for tr in traces_to_remove:
        st.remove(tr)
    if len(st) == 0:
        logger.warning(f"no trace left, station ignored. ({station_id})")
        # print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # attach channel metadata and instrument response
    traces_to_remove = []
    for tr in st:
        try:
            # TODO check if response/metadata changes from tr.stats.starttime to tr.stats.endtime?
            meta = inv.get_channel_metadata(tr.id, tr.stats.starttime)
            tr.stats.metadata = meta
        except Exception as err:
            logger.warning(
                f"cannot get instrument response, {tr.id} ignored. ({station_id})"
            )
            logger.debug(f"{err}")
            # print(
            #     f"{err}\n[WARN] cannot get instrument response, {tr.id} ignored. ({station_id})\n"
            # )
            traces_to_remove.append(tr)
    for tr in traces_to_remove:
        st.remove(tr)
    if len(st) == 0:
        logger.warning(f"no trace left, station ignored. ({station_id})")
        # print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # check if all components have the same coordinates
    for key in ["latitude", "longitude", "elevation", "local_depth"]:
        val0 = st[0].stats.metadata[key]
        if not all([tr.stats.metadata[key] == val0 for tr in st[1:]]):
            logger.warning(
                f"coordinates are not the same across all components, only one is kept. ({station_id})"
            )
            # print(
            #     f"[WARN] coordinates are not the same across all components, only one is kept. ({station_id})\n"
            # )
            break

    # remove trend and apply taper
    st.detrend("linear")
    # st.taper(0.5, max_length=config_taper_width)

    # remove instrument response and resample
    traces_to_remove = []
    t0 = max([tr.stats.starttime for tr in st])
    t1 = min([tr.stats.endtime for tr in st])
    resample_starttime = t0
    resample_npts = int((t1 - t0) * config_sampling_rate)
    for tr in st:
        fs = tr.stats.sampling_rate
        npts = tr.stats.npts
        npad = 5
        nfft = scipy.fft.next_fast_len(npts + 2 * npad)
        freqs = np.fft.rfftfreq(nfft, d=1 / fs)
        # sig_spectrum = np.fft.rfft(tr.data, nfft)
        sig_spectrum = np.fft.rfft(
            np.pad(tr.data, (npad, nfft - npts - npad), "mean"), nfft
        )
        # pre-filter
        lp_sos = scipy.signal.butter(
            config_lowpass_N, config_lowpass_Wn, "lowpass", fs=fs, output="sos"
        )
        # tr.data = scipy.signal.sosfiltfilt(lp_sos, tr.data)
        # tr.data = scipy.signal.sosfiltfilt(hp_sos, tr.data)
        _, h_lp = scipy.signal.freqz_sos(lp_sos, worN=freqs, fs=fs)
        sig_spectrum *= abs(h_lp)
        tr.stats.starttime = tr.stats.starttime - npad / fs
        tr.data = np.fft.irfft(sig_spectrum, nfft)
        # resample by lanczos interpolation
        tr.interpolate(
            config_sampling_rate,
            starttime=resample_starttime,
            npts=resample_npts,
            method="lanczos",
            a=20,
        )
    for tr in traces_to_remove:
        st.remove(tr)
    if len(st) == 0:
        logger.warning(f"no trace left, station ignored. ({station_id})")
        continue

    # check if Z component exists
    if len(st.select(component="Z")) != 1:
        logger.warning(f"no Z component, station ignored. ({station_id})")
        continue

    # check if the dip angle of Z component is -90.0
    tr = st.select(component="Z")[0]
    dip = tr.stats.metadata["dip"]
    if dip != -90.0:  # FIXME what about 90.0?
        logger.warning(
            f"Z component's dip angle ({dip}) != -90, station ignored. ({station_id})"
        )
        # print(
        #     f"[WARN] Z component's dip angle ({dip}) != -90, station ignored. ({station_id})\n"
        # )
        continue

    # check horizontal orientation codes: either 1,2 or E,N
    traces_to_remove = []
    horizontal_orientations = {o for o in orientation_codes if o != "Z"}
    if len(horizontal_orientations) > 0:
        if horizontal_orientations != {"1", "2"} and horizontal_orientations != {
            "E",
            "N",
        }:
            logger.warning(
                f"horizontal components ({horizontal_orientations}) != (1,2) or (N,E), ignored. ({station_id})"
            )
            # print(
            #     f"[WARN] horizontal components ({horizontal_orientations}) != (1,2) or (N,E), ignored. ({station_id})\n"
            # )
            for orientation in horizontal_orientations:
                for tr in st.select(component=orientation):
                    traces_to_remove.append(tr)
        else:
            horizontal_orientations = list(horizontal_orientations)
            st1 = st.select(component=horizontal_orientations[0])
            st2 = st.select(component=horizontal_orientations[1])
            if len(st1) != 1 or len(st2) != 1:
                for tr in st1:
                    traces_to_remove.append(tr)
                for tr in st2:
                    traces_to_remove.append(tr)
            else:
                # check dip/azimuth angles
                tr1, tr2 = st1[0], st2[0]
                dip1, dip2 = tr1.stats.metadata["dip"], tr2.stats.metadata["dip"]
                az1, az2 = tr1.stats.metadata["azimuth"], tr2.stats.metadata["azimuth"]
                if dip1 != 0 or dip2 != 0:
                    logger.warning(
                        f"horizontal components' dip angles ({dip1}, {dip2}) != 0, station ignored. ({station_id})"
                    )
                    # print(
                    #     f"[WARN] horizontal components' dip angles ({dip1}, {dip2}) != 0, ignored. ({station_id})\n"
                    # )
                    traces_to_remove += [tr1, tr2]
                # if abs((az1 - az2)%180 / 90 - 1) > 1e-5:
                if abs((az1 - az2) % 360) < 10:
                    logger.warning(
                        f"horizontal components' azimuth angles ({az1}, {az2}) are close to parallel, station ignored. ({station_id})"
                    )
                    # print(
                    #     f"[WARN] horizontal components' azimuth angles ({az1}, {az2}) are close to parallel, ignored. ({station_id})\n"
                    # )
                    traces_to_remove += [tr1, tr2]
    for tr in traces_to_remove:
        if tr in st:
            st.remove(tr)

    # check if no trace left
    if len(st) == 0:
        logger.warning(f"no trace left, station ignored. ({station_id})")
        # print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    if config["enforce_3comp"] and len(st) != 3:
        logger.warning(f"not exactly 3 components, station ignored. ({station_id})")
        # print(f"[WARN] not exactly 3 components, station ignored. ({station_id})\n")
        continue

    station_name = f"{network}_{station}"
    if station_name in groot:
        logger.warning(f"{station_name} alreay exists, station ignored. ({station_id})")
        # print(f"[WARN] {station_name} alreay exists, overwrite. ({station_id})\n")
        h5f.remove_node(groot, station_name)
    gsta = h5f.create_group(groot, station_name)
    gsta._v_attrs["network"] = st[0].stats.network
    gsta._v_attrs["station"] = st[0].stats.station
    gsta._v_attrs["location"] = st[0].stats.location
    # gsta._v_attrs['band_instrument'] = st[0].stats.channel[0:2]
    gsta._v_attrs["latitude"] = float(st[0].stats.metadata["latitude"])
    gsta._v_attrs["longitude"] = float(st[0].stats.metadata["longitude"])
    gsta._v_attrs["elevation"] = float(st[0].stats.metadata["elevation"])
    gsta._v_attrs["depth"] = float(st[0].stats.metadata["local_depth"])
    gsta._v_attrs["first_arrival"] = first_arrival_time

    # #DEBUG
    # print(f"{resample_starttime._ns=} ")
    # for tr in st:
    #     print(f"{tr.stats.starttime._ns=} ")
    # FIXME trace.interpolate internally represent time as float(ns/1e9) which results in a loss of precision
    #       so the following assertion may fail
    # assert all([tr.stats.starttime == resample_starttime for tr in st])
    new_starttime = st[0].stats.starttime
    assert all([tr.stats.starttime == new_starttime for tr in st])
    assert all([tr.stats.npts == resample_npts for tr in st])
    assert all([tr.stats.sampling_rate == config_sampling_rate for tr in st])

    # store waveforms in array, e.g. /NET_STA/DATA_DISP[0:nchan, 0:npts]
    atom = tables.Atom.from_dtype(np.dtype(np.float32))
    filters = tables.Filters(complevel=3, complib="zlib")
    array_name = "DATA_RAW"
    nchan = len(st)
    shape = (nchan, resample_npts)
    ca = h5f.create_carray(gsta, array_name, atom, shape, filters=filters)
    for i in range(nchan):
        ca[i, :] = np.array(st[i].data, dtype=np.float32)
    # set attributes
    # ca.attrs['network'] = st[0].stats.network
    # ca.attrs['station'] = st[0].stats.station
    # ca.attrs['location'] = st[0].stats.location
    # ca.attrs['latitude'] = float(st[0].stats.metadata['latitude'])
    # ca.attrs['longitude'] = float(st[0].stats.metadata['longitude'])
    # ca.attrs['elevation'] = float(st[0].stats.metadata['elevation'])
    # ca.attrs['depth'] = float(st[0].stats.metadata['local_depth'])
    dtype = [("name", "S3"), ("azimuth", float), ("dip", float)]
    channels = [
        (tr.stats.channel, tr.stats.metadata["azimuth"], tr.stats.metadata["dip"])
        for tr in st
    ]
    ca.attrs["channels"] = np.array(channels, dtype=dtype)
    ca.attrs["starttime"] = (
        new_starttime  # st[0].stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
    )
    ca.attrs["sampling_rate"] = config_sampling_rate
    ca.attrs["npts"] = resample_npts
    ca.attrs["filter"] = [
        {"type": "lowpass", "N": config_lowpass_N, "Wn": config_lowpass_Wn},
    ]
    ca.attrs["type"] = "RAW"
    ca.attrs["unit"] = "count"

    inv_used.extend(inv)
    h5f.flush()

h5f.close()
inv_used.write(out_channel_file, format="stationtxt")

# print("==================================\n")
# print(f"[INFO] END: {datetime.datetime.now()}")
logger.info("END <<<<<<<<<<<<<<<<<<<<<")
