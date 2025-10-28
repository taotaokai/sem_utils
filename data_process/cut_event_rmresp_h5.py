import sys
import os
import datetime
import numpy as np
import tables as tb
import scipy
import yaml

from obspy import read_inventory, Inventory, read_events, geodetics
from obspy.clients.filesystem.sds import Client
from obspy.taup import TauPyModel


""" data structure
/
|- attrs: event, origin_time, location
|
|- NET_STA[0:nchan, 0:npts] : carray
   attrs: net,sta,loc,
          channels[0:nchan], (e.g. [(name:HHE, azimuth:90, dip:0), ('HHN',...), ('HHZ',...)])
          stla,stlo,stel,stdp
          starttime, sampling_rate, npts,
          unit (e.g. count, meter, m/s)
          filter_type, filter_param
          response_corner_frequency

e.g.

/
|- attrs: event_id="C202006211907A", origin_time='2020-06-21T19:07:00.000000', ...
|
|- BE_MEM[0:3, 0:npts]
         |- attrs: network='BE', station='MEM', location='00',
         |         channels=[('HHE',90,0), ('HHN',0,0), ('HHZ',0,-90)],
         |         starttime='2020-03-21T00:49:51.008393',
         |         stla,stlo,stel,stdp
         |         sampling_rate=10.0,
         |         npts=18000,
         |         type='displacement'
         |         unit='meter',
         |         filter_type='butter',
         |         filter_param={'N':8, 'Wn':[0.008, 1]}
         |         response_corner_frequency=[0.01, 5]

NOTE: time precision is microsecond

"""

print("==================================\n")
print(f"[INFO] START: {datetime.datetime.now()}\n")

_DEBUG = True
_DEBUG = False


def _get_resp_corner_freq(freqs, resp_amp, ref_freq=None, threshold_dB=3):
    """get lower/higher corner frequencies of a response curve with flat plateau
    e.g. velocity response for a broadband velocity seismometer
    """
    # ref_idx = np.argmax(amp)
    # ref_amp = amp[ref_idx]
    # amp /= ref_amp # normalize amplitude
    if ref_freq is not None:
        ref_idx = int(np.interp(ref_freq, freqs, np.arange(len(freqs))))
    else:
        ref_idx = np.argmax(resp_amp)
    ref_amp = resp_amp[ref_idx]
    normalized_amp = resp_amp / ref_amp  # normalize amplitude

    assert threshold_dB > 0
    threshold_ratio = 10 ** (-threshold_dB / 20)

    # initial corner frequency
    left_corner_freq, right_corner_freq = 0, np.max(freqs)

    inds = np.flatnonzero(normalized_amp[:ref_idx] <= threshold_ratio)
    if inds.size > 0:
        left_corner_idx = np.max(inds)
        left_corner_freq = freqs[left_corner_idx]
        # left_corner_amp = normalized_amp[left_corner_idx]

    inds = np.flatnonzero(normalized_amp[ref_idx:] <= threshold_ratio)
    if inds.size > 0:
        right_corner_idx = ref_idx + np.min(inds)
        right_corner_freq = freqs[right_corner_idx]

    # if _DEBUG:
    #     print(f'[DEBUG] ref_freq = {freqs[ref_idx]}')
    #     plt.loglog(freqs, normalized_amp)
    #     plt.loglog([ref_freq, ref_freq], [np.min(normalized_amp), np.max(normalized_amp)], 'k')
    #     plt.loglog([left_corner_freq, left_corner_freq], [np.min(normalized_amp), np.max(normalized_amp)], 'r')
    #     plt.loglog([right_corner_freq, right_corner_freq], [np.min(normalized_amp), np.max(normalized_amp)], 'b')
    #     plt.show()

    return left_corner_freq, right_corner_freq


#
config_yaml = sys.argv[1]
CMT_file = sys.argv[2]
out_h5_file = sys.argv[3]
out_channel_file = sys.argv[4]

if os.path.dirname(out_h5_file):
    os.makedirs(os.path.dirname(out_h5_file), exist_ok=True)
if os.path.dirname(out_channel_file):
    os.makedirs(os.path.dirname(out_channel_file), exist_ok=True)

# h5 file handle
h5f = tb.open_file(out_h5_file, mode="w")

# event
event_data = read_events(CMT_file)[0]
event_id = event_data.event_descriptions[0].text
event_origin = event_data.preferred_origin()

print("==================================\n")
print(f"[INFO] event: {event_origin}\n")

# config
with open(config_yaml, "r") as file:
    config = yaml.safe_load(file)

config_zcomp_pattern = config["zcomp_pattern"]

config_sampling_rate = config["resample"]["sampling_rate"]
config_lowpass_N = config["resample"]["filter"]["N"]
config_lowpass_Wn = config["resample"]["filter"]["Wn"]

config_output_type = config["remove_response"]["output_type"]
config_output_unit = config["remove_response"]["output_unit"]
config_ref_resp_type = config["remove_response"]["ref_resp_type"]
config_ref_freq = config["remove_response"]["ref_freq"]
config_ref_threshold = config["remove_response"]["ref_threshold"]
config_max_resp_lc = config["remove_response"]["max_resp_lc"]
config_min_resp_hc = config["remove_response"]["min_resp_hc"]
# config_filter_lstop = config['remove_response']['filter']['lstop']
# config_filter_min_lpass = config['remove_response']['filter']['min_lpass']
config_highpass_N = config["remove_response"]["filter"]["N"]
config_highpass_min_cutoff_freq = config["remove_response"]["filter"]["min_cutoff_freq"]
# config_filter_max_hpass = config['remove_response']['filter']['max_hpass']
# config_filter_gpass = config['remove_response']['filter']['gpass']
# config_filter_gstop = config['remove_response']['filter']['gstop']

time_before_first_arrival = config["time_window"]["before_first_arrival_seconds"]
# time_after_first_arrival = config['data']['time_window']['after_first_arrival_seconds']
time_after_origin = config["time_window"]["after_origin_time_seconds"]

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

t0 = event_origin.time - time_before_first_arrival
t1 = event_origin.time + time_after_origin

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
        print(f"[WARN] fail to read {sxml_file}, skip. (Error: {err})")
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
            print(f"[WARN] fail to select inventory for {nslc} on {event_origin.time}")
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
    print("==================================")
    station_id = ".".join([network, station, location, channel[:-1] + "?"])
    print(f"[INFO] {station_id}\n")

    # read StationXML file
    sxml_file = os.path.join(stations_path, f"{network}.{station}.xml")
    try:
        inv = read_inventory(sxml_file)
    except Exception as err:
        print(
            f"{err}\n[WARN] fail to read {sxml_file}, station ignored. ({station_id})\n"
        )
        continue

    inv = inv.select(
        network=network,
        station=station,
        location=location,
        channel=channel[:-1] + "?",
        time=event_origin.time,
    )
    if not inv:
        print(
            f"[WARN] no station metadata on {event_origin.time}), station ignored. ({station_id})\n"
        )
        continue

    # get channel coordinates
    try:
        seed_id = ".".join([network, station, location, channel])
        station_coords = inv.get_coordinates(seed_id, event_origin.time)
    except Exception as err:
        print(
            f"{err}\n[WARN] cannot find station metadata for {seed_id}, station ignored. ({station_id})\n"
        )
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

    # data time window
    time_window_starttime = first_arrival_time - time_before_first_arrival
    time_window_endtime = event_origin.time + time_after_origin
    # round time to integer seconds
    if time_window_starttime.microsecond != 0:
        time_window_starttime = time_window_starttime.replace(microsecond=0)
    if time_window_endtime.microsecond != 0:
        time_window_endtime = (time_window_endtime + 1).replace(microsecond=0)

    # get waveforms
    try:
        pad_time = 5  # seconds
        st = client.get_waveforms(
            network,
            station,
            location,
            channel[0:2] + "?",
            time_window_starttime - pad_time,
            time_window_endtime + pad_time,
            merge=-1,
        )
    except Exception as err:
        print(
            f"{err}\n[WARN] failed to get waveforms, station ignored. ({station_id})\n"
        )
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
            print(f"{st1.__str__(extended=True)}")
            print(f"[WARN] gaps/overlaps found, {st1[0].id} ignored. ({station_id})\n")
            for tr in st1:
                st.remove(tr)
    if len(st) == 0:
        print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # check time coverage
    traces_to_remove = []
    for tr in st:
        dt = tr.stats.delta
        if tr.stats.starttime > time_window_starttime:
            print(
                f"{tr}\n[WARN] trace starttime {tr.stats.starttime} is later than required {time_window_starttime}, {tr.id} ignored. ({station_id})\n"
            )
            traces_to_remove.append(tr)
        elif tr.stats.endtime < time_window_endtime:
            print(
                f"{tr}\n[WARN] trace endtime {tr.stats.endtime} is earlier than required {time_window_endtime}, {tr.id} ignored. ({station_id})\n"
            )
            traces_to_remove.append(tr)
    for tr in traces_to_remove:
        st.remove(tr)
    if len(st) == 0:
        print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # attach channel metadata and instrument response
    traces_to_remove = []
    for tr in st:
        try:
            # TODO check if response/metadata changes from tr.stats.starttime to tr.stats.endtime?
            meta = inv.get_channel_metadata(tr.id, tr.stats.starttime)
            resp = inv.get_response(tr.id, tr.stats.starttime)
            tr.stats.metadata = meta
            tr.stats.response = resp
        except Exception as err:
            print(
                f"{err}\n[WARN] cannot get instrument response, {tr.id} ignored. ({station_id})\n"
            )
            traces_to_remove.append(tr)
    for tr in traces_to_remove:
        st.remove(tr)
    if len(st) == 0:
        print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # check if all components have the same coordinates
    for key in ["latitude", "longitude", "elevation", "local_depth"]:
        val0 = st[0].stats.metadata[key]
        if not all([tr.stats.metadata[key] == val0 for tr in st[1:]]):
            print(
                f"[WARN] coordinates are not the same across all components, only one is kept. ({station_id})\n"
            )
            break

    # get response corner frequency on response spectrum
    traces_to_remove = []
    # frequency samples [1/1000, 1/999, ... , 0.5, 1, 1.1, 1.2, ...] Hz
    ref_freqs = np.hstack(
        (1.0 / np.arange(1000, 1, -1), np.arange(1, config_sampling_rate, 0.1))
    )
    for tr in st:
        fs = tr.stats.sampling_rate
        try:
            # get instrument response spectrum
            resp = tr.stats.response
            ref_resp = resp.get_evalresp_response_for_frequencies(
                ref_freqs, output=config_ref_resp_type
            )
        except Exception as err:
            print(
                f"{err}\n[WARN] failed to evaluate instrument response, {tr.id} ignored. ({station_id})\n"
            )
            traces_to_remove.append(tr)
            continue
        # print(resp.get_paz())
        flc, fhc = _get_resp_corner_freq(
            ref_freqs,
            np.abs(ref_resp),
            ref_freq=config_ref_freq,
            threshold_dB=config_ref_threshold,
        )
        tr.stats.resposne_corner_frequency = [flc, fhc]
        # number of zeros
        paz = resp.get_paz()
        nz = sum([1 if z == 0 else 0 for z in paz.zeros])
        print(f"num. of zero zeros = {nz}")

    for tr in traces_to_remove:
        st.remove(tr)
    if len(st) == 0:
        print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # determine common response frequency band for all components
    resp_lc = max(
        [
            tr.stats.resposne_corner_frequency[0]
            for tr in st
            if "resposne_corner_frequency" in tr.stats
        ]
    )
    resp_hc = min(
        [
            tr.stats.resposne_corner_frequency[1]
            for tr in st
            if "resposne_corner_frequency" in tr.stats
        ]
    )
    print(
        f"[INFO] response corner frequency = [{resp_lc}, {resp_hc}] for {station_id}\n"
    )
    if resp_lc > config_max_resp_lc:
        print(
            f"[WARN] response lower corner frequency > {config_max_resp_lc} Hz, station ignored. ({station_id})\n"
        )
        continue
    if resp_hc < config_min_resp_hc:
        print(
            f"[WARN] response higher corner frequency < {config_min_resp_hc} Hz, station ignored. ({station_id})\n"
        )
        continue

    # # butter-worth filter pass/stop band
    # lstop = config_filter_lstop
    # lpass = max(resp_lc, config_filter_min_lpass)
    # hpass = min(resp_hc, config_filter_max_hpass)
    # nyq_freq = config_sampling_rate * 0.5
    # hstop = nyq_freq
    # if lpass >= hpass:
    #     print(f'[WARN] lpass({lpass}) >= hpass({hpass}), station ignored. ({station_id})\n')
    #     continue
    # max_fs = max([tr.stats.sampling_rate for tr in st])
    # butter_N, butter_Wn = scipy.signal.buttord([lpass, hpass], [lstop, hstop], config_filter_gpass, config_filter_gstop, fs=max_fs)
    # print(f'[INFO] filter design pass/stop band: [{lpass}, {hpass}], [{lstop}, {hstop}]\n')
    # print(f'[INFO] filter butter: N, Wn = {butter_N}, {butter_Wn}\n')
    config_highpass_Wn = max(config_highpass_min_cutoff_freq, resp_lc)

    # remove trend and apply taper
    st.detrend("linear")
    # st.taper(0.5, max_length=config_taper_width)

    # for debug
    if _DEBUG:
        st1 = st.copy()
        for tr in st1:
            tr.stats.network = "RAW_" + tr.stats.network

    # remove instrument response and resample
    traces_to_remove = []
    resample_starttime = time_window_starttime
    resample_npts = int(
        (time_window_endtime - time_window_starttime) * config_sampling_rate
    )
    for tr in st:
        fs = tr.stats.sampling_rate
        npts = tr.stats.npts
        npad = int(2.0 * fs / config_highpass_Wn)
        # npad = int(2.0 / min(config_lower_corner_taper_width, config_higher_corner_taper_width) * fs)
        # print(f'[INFO] resample npad = {npad}')
        nfft = scipy.fft.next_fast_len(npts + npad)
        freqs = np.fft.rfftfreq(nfft, d=1 / fs)
        sig_spectrum = np.fft.rfft(tr.data, nfft)
        # pre-filter
        lp_sos = scipy.signal.butter(
            config_lowpass_N, config_lowpass_Wn, "lowpass", fs=fs, output="sos"
        )
        hp_sos = scipy.signal.butter(
            config_highpass_N, config_highpass_Wn, "highpass", fs=fs, output="sos"
        )
        _, h_lp = scipy.signal.freqz_sos(lp_sos, worN=freqs, fs=fs)
        _, h_hp = scipy.signal.freqz_sos(hp_sos, worN=freqs, fs=fs)
        sig_spectrum *= abs(h_lp) * abs(h_hp)
        # if _DEBUG:
        #     plt.loglog(freqs, abs(h_lp) * abs(h_hp))
        #     plt.title('Butterworth filter frequency response')
        #     plt.xlabel('Frequency [Hz]')
        #     plt.ylabel('Amplitude')
        #     plt.grid(which='both', axis='both')
        #     plt.show()
        # remove instrument response
        try:
            resp = tr.stats.response
            output_resp = resp.get_evalresp_response_for_frequencies(
                freqs, output=config_output_type
            )
        except Exception as err:
            print(
                f"{err}\n[WARN] failed to evaluate instrument response, {tr.id} ignored. ({station_id})\n"
            )
            traces_to_remove.append(tr)
            continue
        inds = (freqs != 0) & (abs(output_resp) > 0)
        sig_spectrum[inds] /= output_resp[inds]
        sig_spectrum[0] = 0  # ensure no amplitude at zero-frequency
        tr.data = np.fft.irfft(sig_spectrum, nfft)[:npts]
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
        print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    # st.taper(0.5, max_length=config_taper_width)

    # for debug
    if _DEBUG:
        st1.extend(st)
        # st1.filter('lowpass', freq=0.02)
        st1.plot(equal_scale=False)

    # check if Z component exists
    if len(st.select(component="Z")) != 1:
        print(f"[WARN] no Z component, station ignored. ({station_id})\n")
        continue

    # check if the dip angle of Z component is -90.0
    tr = st.select(component="Z")[0]
    dip = tr.stats.metadata["dip"]
    if dip != -90.0:  # FIXME what about 90.0?
        print(
            f"[WARN] Z component's dip angle ({dip}) != -90, station ignored. ({station_id})\n"
        )
        continue

    # check horizontal orientation codes: either 1,2 or E,N
    traces_to_remove = []
    horizontal_orientations = {o for o in orientation_codes if o != "Z"}
    if len(horizontal_orientations) > 0:
        if horizontal_orientations != {"1", "2"} and horizontal_orientations != {
            "E",
            "N",
        }:
            print(
                f"[WARN] horizontal components ({horizontal_orientations}) != (1,2) or (N,E), ignored. ({station_id})\n"
            )
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
                    print(
                        f"[WARN] horizontal components' dip angles ({dip1}, {dip2}) != 0, ignored. ({station_id})\n"
                    )
                    traces_to_remove += [tr1, tr2]
                # if abs((az1 - az2)%180 / 90 - 1) > 1e-5:
                if abs((az1 - az2) % 360) < 10:
                    print(
                        f"[WARN] horizontal components' azimuth angles ({az1}, {az2}) are close to parallel, ignored. ({station_id})\n"
                    )
                    traces_to_remove += [tr1, tr2]
    for tr in traces_to_remove:
        if tr in st:
            st.remove(tr)

    # check if no trace left
    if len(st) == 0:
        print(f"[WARN] no trace left, station ignored. ({station_id})\n")
        continue

    if config["enforce_3comp"] and len(st) != 3:
        print(f"[WARN] not exactly 3 components, station ignored. ({station_id})\n")
        continue

    station_name = f"{network}_{station}"
    if station_name in groot:
        print(f"[WARN] {station_name} alreay exists, overwrite. ({station_id})\n")
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

    assert all([tr.stats.starttime == resample_starttime for tr in st])
    assert all([tr.stats.npts == resample_npts for tr in st])
    assert all([tr.stats.sampling_rate == config_sampling_rate for tr in st])

    # store waveforms in array, e.g. /NET_STA/DATA_DISP[0:nchan, 0:npts]
    atom = tb.Atom.from_dtype(np.dtype(np.float32))
    filters = tb.Filters(complevel=3, complib="zlib")
    array_name = f"DATA_{config_output_type}"
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
        resample_starttime  # st[0].stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
    )
    ca.attrs["sampling_rate"] = config_sampling_rate
    ca.attrs["npts"] = resample_npts
    ca.attrs["filter"] = [
        {"type": "lowpass", "N": config_lowpass_N, "Wn": config_lowpass_Wn},
        {"type": "highpass", "N": config_highpass_N, "Wn": config_highpass_Wn},
    ]
    ca.attrs["response_corner_frequency"] = [resp_lc, resp_hc]
    ca.attrs["type"] = config_output_type
    ca.attrs["unit"] = config_output_unit

    inv_used.extend(inv)
    h5f.flush()

h5f.close()
inv_used.write(out_channel_file, format="stationtxt")

print("==================================\n")
print(f"[INFO] END: {datetime.datetime.now()}")
