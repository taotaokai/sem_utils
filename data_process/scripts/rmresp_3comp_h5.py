import sys
import os
import glob
import datetime
import numpy as np
import tables as pt
import scipy
import matplotlib.pyplot as plt
import yaml

from obspy.signal.invsim import cosine_sac_taper
from obspy import read, read_inventory, Inventory, read_events, geodetics
# from obspy.signal.util import _npts2nfft
from obspy.clients.filesystem.tsindex import Client
from obspy.taup import TauPyModel

""" data structure
/STATION_<net>_<sta>/TRACE_<starttime>_<endtime>[0:ncomp, 0:npts]
                     attrs: orientation[0:ncomp]: [(code, azimuth, dip), ...]
                            starttime, sampling_rate, npts,
                            data_type(e.g. RAW, DISP, VEL),
                            unit(e.g. count, meter, m/s)
                            filter_type, filter_limit
                            response_cutoff_frequency

e.g.

/STATION_BE_MEM/TRACE_20200321T004951_20200321T011950[3,npts]
                attrs: starttime='2020-03-21T00:49:51.008393',
                       sampling_rate=10.0,
                       npts=18000,
                       orientation=[('E', 90, 0),('N', 0, 0),('Z', 0, -90)]
                       data_type='DISP',
                       data_unit='meter',
                       filter_limit=[0.01, 0.02, 10, 12],
                       filter_type='cosine_sac_taper'
                       response_cutoff_frequency=[0.01, 5]

NOTE: time precision is microsecond

"""

print('==================================\n')
print(f'[INFO] START: {datetime.datetime.now()}\n')

_DEBUG = True
_DEBUG = False

def _get_cutoff_freq(freqs, amp, ref_freq=-1, cutoff_threshold=0.05):
    """ get lower/higher cutoff frequencies of a response curve with flat top
        e.g. velocity response for a broadband velocity seismometer
    """
    # ref_idx = np.argmax(amp)
    # ref_amp = amp[ref_idx]
    # amp /= ref_amp # normalize amplitude
    if ref_freq > 0:
        ref_idx = int(np.interp(ref_freq, freqs, np.arange(len(freqs))))
    else:
        ref_idx = np.argmax(amp)
    ref_amp = amp[ref_idx]
    normalized_amp = amp / ref_amp # normalize amplitude

    # initial cutoff frequency
    left_cutoff_freq, right_cutoff_freq = 0, np.max(freqs)

    inds = np.flatnonzero(normalized_amp[:ref_idx] <= cutoff_threshold)
    if inds.size > 0:
        left_cutoff_idx = np.max(inds)
        left_cutoff_freq = freqs[left_cutoff_idx]
        # left_cutoff_amp = normalized_amp[left_cutoff_idx]

    inds = np.flatnonzero(normalized_amp[ref_idx:] <= cutoff_threshold)
    if inds.size > 0:
        right_cutoff_idx = ref_idx + np.min(inds)
        right_cutoff_freq = freqs[right_cutoff_idx]

    if _DEBUG:
        print(f'[DEBUG] ref_freq = {freqs[ref_idx]}')
        plt.loglog(freqs, normalized_amp)
        plt.loglog([ref_freq, ref_freq], [np.min(normalized_amp), np.max(normalized_amp)], 'k')
        plt.loglog([left_cutoff_freq, left_cutoff_freq], [np.min(normalized_amp), np.max(normalized_amp)], 'r')
        plt.loglog([right_cutoff_freq, right_cutoff_freq], [np.min(normalized_amp), np.max(normalized_amp)], 'b')
        plt.show()

    return left_cutoff_freq, right_cutoff_freq

# config_yaml = 'config.yaml'
# CMT_file = 'CMTSOLUTION'
# h5_file = 'data_remove_resp.h5'

config_yaml = sys.argv[1]
CMT_file = sys.argv[2]
h5_file = sys.argv[3]

# h5 handle
h5f = pt.open_file(h5_file, mode="w")

# event
event_data = read_events(CMT_file)[0]
event_origin = event_data.preferred_origin()

print('==================================\n')
print(f'[INFO] event: {event_origin}\n')

# config
with open(config_yaml, 'r') as file:
    config = yaml.safe_load(file)

sampling_rate = config['data']['sampling_rate_Hz']

resp_ref_freq = config['data']['remove_response']['ref_freq_Hz']
resp_cutoff_threshold = config['data']['remove_response']['cutoff_threshold']
max_left_cutoff_frequency = config['data']['remove_response']['max_left_cutoff_frequency_Hz']
min_right_cutoff_frequency = config['data']['remove_response']['min_right_cutoff_frequency_Hz']

# lc_taper, rc_taper: frequency domian cosine window [flc, flc + lc_taper, frc - rc_taper, frc]
#                     will introduce side lope of duration ~ 1/taper seconds on both sides of the spike response
lc_taper = config['data']['remove_response']['left_cutoff_taper_width_Hz']
rc_taper = config['data']['remove_response']['right_cutoff_taper_width_Hz']

time_before_first_arrival = config['data']['time_window']['before_first_arrival_seconds']
# time_after_first_arrival = config['data']['time_window']['after_first_arrival_seconds']
time_after_origin = config['data']['time_window']['after_origin_time_seconds']

# for merging traces
misalignment_threshold = config['data']['misalignment_threshold']

stations_path = config['data']['stations_path']
waveforms_path = config['data']['waveforms_path']
indexdb_path = config['data']['indexdb_path']

client = Client(indexdb_path, datapath_replace=("^", os.path.join(waveforms_path,'')))

# read in all staitonxml files
inv = Inventory()
xml_list = glob.glob(os.path.join(stations_path, '*.xml'))
for xml_file in xml_list:
    try:
        inv1 = read_inventory(xml_file)
    except Exception as err:
        print(f'{err}\n[WARN] failed to read {xml_file}\n')
        continue
    inv.extend(inv1)

# taup
taup_model = TauPyModel(model="ak135")

# query for Z component
extents_zcomp = client.get_availability_extent(channel="?HZ")
# DEBUG
# net, sta = 'DK', 'STE05'
# net, sta = 'NO', 'BRBA'
# net, sta = 'NO', 'NC405'
# net, sta = 'C4', 'CERNS'
# extents_zcomp = client.get_availability_extent(network=net, station=sta, channel="?HZ")

# keep only one 3-comp data for each (net, sta)
# { (net, sta): [(loc, band, extent), ...], ... }
print('==================================\n')
net_sta_dict = {}
for extent in extents_zcomp:
    network = extent[0]
    station = extent[1]
    location = extent[2]
    channel = extent[3]
    band = channel[0]
    net_sta = (network, station)
    loc_band = (location, band, extent)
    if net_sta not in net_sta_dict:
        net_sta_dict[net_sta] = [loc_band]
    else:
        net_sta_dict[net_sta].append(loc_band)

band_code_priority = config['data']['band_code_priority']
loc_code_priority = config['data']['location_code_priority']
extents_used = []
for net_sta in net_sta_dict:
    loc_band_list = net_sta_dict[net_sta]
    # get index of band/loc code according to order in band/loc_code_order
    n = max(len(band_code_priority), len(loc_code_priority))
    inds_band = [band_code_priority.index(x[1]) if x[1] in band_code_priority else n for x in loc_band_list]
    inds_loc = [loc_code_priority.index(x[0]) if x[0] in loc_code_priority else n for x in loc_band_list]
    sort_idx = np.lexsort((inds_loc, inds_band))
    if '*' in loc_code_priority:
        inds_loc = [x if x < n else -1 for x in inds_loc]
    flag_found = False
    for idx in sort_idx:
        if inds_band[idx] != n and inds_loc[idx] != n:
            flag_found = True
            extents_used.append(loc_band_list[idx][2])
            break
    if not flag_found:
        print(f'[WARN] skip {net_sta}: {loc_band_list}')
        continue

del net_sta_dict

print(f'[INFO] stations to process: \n')
for ext in extents_used:
    print(f'{ext}')
print(f'\n')

# loop all stations
# for network, station, location, channel, begintime, endtime in extents_zcomp:
for network, station, location, channel, begintime, endtime in extents_used:
    print('==================================')
    station_id = '.'.join([network, station, location, channel[:-1]+'?'])
    print(f'[INFO] {station_id}\n')

    # get channel coordinates
    try:
        seed_id = '.'.join([network, station, location, channel])
        station_coords = inv.get_coordinates(seed_id, event_origin.time)
    except Exception as err:
        print(f"{err}\n[WARN] cannot find station metadata for {seed_id}, ignore station {station_id}\n")
        continue

    # calculate event-station distance
    dist, az, baz = geodetics.gps2dist_azimuth(
        event_origin.latitude, event_origin.longitude,
        station_coords['latitude'], station_coords['longitude'])
    dist_degree = geodetics.kilometer2degrees(dist/1000.0)

    # calculate first arrival time using ak135 model
    arrivals = taup_model.get_travel_times(
        source_depth_in_km=event_origin.depth/1000.0,
        distance_in_degree=dist_degree)
    first_arrival_time = event_origin.time + min([arr.time for arr in arrivals])
    # print(dist, dist_degree, arrivals)

    # data time window
    time_window_starttime = first_arrival_time - time_before_first_arrival
    # time_window_endtime = first_arrival_time + time_after_first_arrival
    time_window_endtime = event_origin.time + time_after_origin
    # print(first_arrival_time, time_window_starttime, time_window_endtime)

    # get waveforms
    try:
        st = client.get_waveforms(network, station, location, channel[0:2]+'?', time_window_starttime, time_window_endtime, merge=False)
    except Exception as err:
        print(f"{err}\n[WARN] failed to get waveforms, ignore station {station_id}\n")
        continue

    # merge adjacent traces with misalignment less than sub-sampling interval
    st.merge(-1, misalignment_threshold=misalignment_threshold)

    # ignore traces with gaps or overlaps
    orientation_codes = set([tr.stats.channel[-1] for tr in st])
    for orientation in orientation_codes:
        st1 = st.select(component=orientation)
        if len(st1) != 1:
            print(f'{st1.__str__(extended=True)}\n[WARN] gaps/overlaps found, traces {st1[0].id} are removed\n')
            for tr in st1: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace is left, ignore station {station_id}\n')
        continue

    # check time coverage
    traces_to_remove = []
    for tr in st:
        dt = tr.stats.delta
        if tr.stats.starttime - time_window_starttime > dt:
            print(f'{tr}\n[WARN] trace\'s starttime {tr.stats.starttime} is later than required {time_window_starttime} by over one sampling interval, trace is ignored\n')
            traces_to_remove.append(tr)
        elif tr.stats.endtime - time_window_endtime > dt:
            print(f'{tr}\n[WARN] trace\'s endtime {tr.stats.endtime} is earlier than required {time_window_endtime} by over one sampling interval, trace is ignored\n')
            traces_to_remove.append(tr)
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace is left, ignore station {station_id}\n')
        continue

#     # check if Z component exists
#     orientation_codes = [tr.stats.channel[-1] for tr in st]
#     if 'Z' not in orientation_codes:
#         print(f'[WARN] Z component not found, station {station_id} is ignored\n')
#         continue
#
#     # check horizontal orientation codes: either 1,2 or E,N
#     horizontal_orientations = {o for o in orientation_codes if o != 'Z'}
#     if horizontal_orientations != {'1', '2'} and horizontal_orientations != {'E', 'N'}:
#         print(f'[WARN] unexpected horizontal orientations {horizontal_orientations}, should be either 1,2 or N,E, horizontal components ignored\n')
#         for orientation in horizontal_orientations:
#             for tr in st.select(component=orientation):
#                 st.remove(tr)

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
            print(f"{err}\n[WARN] cannot get instrument response for {tr.id}, skipped.\n")
            traces_to_remove.append(tr)
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace is left, ignore station {station_id}\n')
        continue

    # check if all components have the same coordinates
    for key in ['latitude','longitude','elevation','local_depth']:
        val0 = st[0].stats.metadata[key]
        if not all([tr.stats.metadata[key] == val0 for tr in st[1:]]):
            print(f'[WARN] coordinates are not the same across all components\n')
            break

    # remove trend
    st.detrend('linear')

    # determine common time range for all components before resample
    resample_starttime = max([tr.stats.starttime for tr in st])
    endtime = min([tr.stats.endtime for tr in st])
    if resample_starttime.microsecond != 0:
        resample_starttime = (resample_starttime + 1).replace(microsecond=0) # round to integer seconds
    resample_npts = int((endtime - resample_starttime) * sampling_rate)

    # resample
    nyq_freq = sampling_rate * 0.5
    resample_hstop = nyq_freq
    resample_hpass = resample_hstop - rc_taper
    resample_prefilt = [-2, -1, resample_hpass, resample_hstop]
    print(f'[INFO] resample_prefilt = {resample_prefilt}\n')
    for tr in st:
        # apply lowpass filter before resampling
        fs = tr.stats.sampling_rate
        npts = tr.stats.npts
        npad = int(2.0 / rc_taper * fs)
        # print(f'[INFO] resample npad = {npad}')
        nfft = scipy.fft.next_fast_len(npts + npad)
        freqs = np.fft.rfftfreq(nfft, d=1/fs)
        #
        sig_spectrum = np.fft.rfft(tr.data, nfft)
        # print(f'[INFO] nyquist frequency = {nyq_freq}')
        # print(f'[INFO] resample_prefilt = {resample_prefilt}')
        taper = cosine_sac_taper(freqs, resample_prefilt) # apply taper near the end
        sig_spectrum *= taper
        tr.data = np.fft.irfft(sig_spectrum, nfft)[:npts]
        # resample by lanczos interpolation
        tr.interpolate(sampling_rate, starttime=resample_starttime, npts=resample_npts, method='lanczos', a=20)

    # for debug
    if _DEBUG:
        st1 = st.copy()
        for tr in st1:
            tr.stats.network = 'RAW_' + tr.stats.network

    # get response corner frequency on velocity response spectrum
    traces_to_remove = []
    for tr in st:
        fs = tr.stats.sampling_rate
        # frequency samples [1/1000, 1/999, ... , 0.5, 1, 1.5, ...] Hz
        raw_freqs = np.hstack((1.0/np.arange(1000, 1, -1), np.arange(1, fs/2+rc_taper+1, 0.5)))
        try:
            # get response for velocity
            resp = tr.stats.response
            raw_resp = resp.get_evalresp_response_for_frequencies(raw_freqs, output='VEL')
        except Exception as err:
            print(f"{err}\n[WARN] failed to evaluate instrument response for {tr.id}, skipped.\n")
            traces_to_remove.append(tr)
            continue
        # # get corner frequency from velocity response spectrum
        # try:
        #     # ref_freq = resp.response_stages[0].normalization_frequency
        #     ref_freq = resp.response_stages[0].stage_gain_frequency
        # except Exception as err:
        #     print(f"{err}\n[WARN] cannot get response normalization frequency for {tr.id}, set to 1.0 Hz.\n")
        #     ref_freq = 1.0
        # print(f'[INFO] ref_freq = {ref_freq}')
        flc, frc = _get_cutoff_freq(raw_freqs, np.abs(raw_resp), ref_freq=resp_ref_freq, cutoff_threshold=resp_cutoff_threshold)
        # flc, frc = _get_cutoff_freq(raw_freqs, np.abs(raw_resp), cutoff_threshold=resp_cutoff_threshold)
        # print(f'[INFO] response cutoff frequency = [{flc}, {frc}] for {tr.id}')
        tr.stats.resposne_corner_frequency = [flc, frc]
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace is left, ignore station {station_id}\n')
        continue

    # determine common frequency range used in removing response for all components
    max_flc = max([tr.stats.resposne_corner_frequency[0] for tr in st if 'resposne_corner_frequency' in tr.stats])
    min_frc = min([tr.stats.resposne_corner_frequency[1] for tr in st if 'resposne_corner_frequency' in tr.stats])
    print(f'[INFO] response cutoff frequency = [{max_flc}, {min_frc}] for {station_id}\n')
    if max_flc > max_left_cutoff_frequency:
        print(f'[WARN] response left cutoff frequency is larger than required {max_left_cutoff_frequency}, ignore station {station_id}\n')
        continue
    if min_frc < min_right_cutoff_frequency:
        print(f'[WARN] response right cutoff frequency is smaller than required {min_right_cutoff_frequency}, ignore station {station_id}\n')
        continue
    rmresp_lstop = max_flc
    rmresp_lpass = max_flc + lc_taper # filter half duration about 1.0/lc_taper sec
    if max_flc <= 0: # no cutoff on low frequency
        rmresp_lstop = -2
        rmresp_lpass = -1
    min_rstop = min(resample_hpass, min_frc)
    rmresp_rstop = min_rstop
    rmresp_rpass = min_rstop - rc_taper
    if rmresp_lpass >= rmresp_rpass:
        print(f'[WARN] rmresp_lpass({rmresp_lpass}) >= rmresp_rpass({rmresp_rpass}), ignore station {station_id}\n')
        continue
    rmresp_prefilt = [rmresp_lstop, rmresp_lpass, rmresp_rpass, rmresp_rstop]
    print(f'[INFO] rmresp_prefilt = {rmresp_prefilt}\n')

    # remove response
    traces_to_remove = []
    for tr in st:
        # apply taper in frequency domain
        npts = tr.stats.npts
        npad = int(2.0 / min(lc_taper, rc_taper) * fs) # approximate length of the bandpass filter
        # print(f'[INFO] rmresp npad = {npad}')
        nfft = scipy.fft.next_fast_len(npts + npad)
        # nfft = _npts2nfft(npts)
        freqs = np.fft.rfftfreq(nfft, d=1/fs)
        sig_spectrum = np.fft.rfft(tr.data, nfft)
        taper = cosine_sac_taper(freqs, rmresp_prefilt)
        sig_spectrum *= taper
        try:
            resp = tr.stats.response
            disp_resp = resp.get_evalresp_response_for_frequencies(freqs, output='DISP')
        except Exception as err:
            print(f"{err}\n[WARN] failed to evaluate instrument response for {tr.id}, trace is skipped.\n")
            traces_to_remove.append(tr)
            continue
        inds = (freqs != 0) & (freqs >= rmresp_lstop) & (freqs <= rmresp_rstop)
        sig_spectrum[inds] /= disp_resp[inds]
        sig_rmresp = np.fft.irfft(sig_spectrum, nfft)[:npts]
        tr.data = sig_rmresp
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace is left, ignore station {station_id}\n')
        continue

    # for debug
    if _DEBUG:
        (st1 + st).plot(equal_scale=False)

    # check if Z component exists
    if len(st.select(component='Z')) != 1:
        print(f'[WARN] Z component not found, ignore station {station_id}\n')
        continue

    # check if the dip angle of Z component is -90.0
    tr = st.select(component='Z')[0]
    dip = tr.stats.metadata['dip']
    if dip != -90.0:
        print(f'[WARN] dip angle of the Z component is {dip}, should be -90! ignore station {station_id}\n')
        continue

    # check horizontal orientation codes: either 1,2 or E,N
    traces_to_remove = []
    horizontal_orientations = {o for o in orientation_codes if o != 'Z'}
    if horizontal_orientations != {'1', '2'} and horizontal_orientations != {'E', 'N'}:
        print(f'[WARN] unexpected horizontal orientations {horizontal_orientations}, should be either 1,2 or N,E, ignore horizontal components\n')
        for orientation in horizontal_orientations:
            for tr in st.select(component=orientation):
                traces_to_remove.append(tr)
    else:
        horizontal_orientations = list(horizontal_orientations)
        st1 = st.select(component=horizontal_orientations[0])
        st2 = st.select(component=horizontal_orientations[1])
        if len(st1) != 1 or len(st2) != 1:
            for tr in st1: traces_to_remove.append(tr)
            for tr in st2: traces_to_remove.append(tr)
        else:
            # check dip/azimuth angles
            tr1, tr2 = st1[0], st2[0]
            dip1, dip2 = tr1.stats.metadata['dip'], tr2.stats.metadata['dip']
            az1, az2 = tr1.stats.metadata['azimuth'], tr2.stats.metadata['azimuth']
            if dip1 != 0 or dip2 != 0:
                print(f'[WARN] {station_id}: horizontal components have dips: {dip1}, {dip2}, should be 0, ignore horizontal components\n')
                traces_to_remove += [tr1, tr2]
            if (az1 - az2) % 180 != 90:
                print(f'[WARN] {station_id}: horizontal components have azimuths: {az1}, {az2}, should be orthogonal, ignore horizontal components\n')
                traces_to_remove += [tr1, tr2]
    for tr in traces_to_remove: st.remove(tr)

    # check if no trace left
    if len(st) == 0:
        print(f'[WARN] no trace is left, ignore station {station_id}\n')
        continue

    if config['data']['enforce_3comp'] and len(st) != 3:
        print(f'[WARN] not exactly 3 components, ignore station {station_id}\n')
        continue

    # create group if not exists: e.g. /STATION_BE_MEM/
    grp_name = '_'.join([network, station]) # net,sta
    grp_name = 'STATION_' + grp_name
    if grp_name not in h5f.root:
        grp = h5f.create_group(h5f.root, grp_name)
    else:
        print(f'[WARN] {grp_name} already exists, overwrite attributes!\n')
        grp = h5f.get_node(h5f.root, grp_name)

    # overwrite attributes
    grp._v_attrs['network'] = network
    grp._v_attrs['station'] = station
    grp._v_attrs['location'] = location
    grp._v_attrs['band'] = channel[0]
    grp._v_attrs['instrument'] = channel[1]
    grp._v_attrs['latitude'] = st[0].stats.metadata['latitude']
    grp._v_attrs['longitude'] = st[0].stats.metadata['longitude']
    grp._v_attrs['elevation'] = st[0].stats.metadata['elevation']
    grp._v_attrs['local_depth'] = st[0].stats.metadata['local_depth']

    # store waveforms in array, e.g. /BE_MEM_00_H_H/DISP_20200321T004951_20200321T011950[3,npts]
    atom = pt.Atom.from_dtype(np.dtype(np.float32))
    filters = pt.Filters(complevel=5, complib='zlib')
    starttime = st[0].stats.starttime.strftime('%Y%m%dT%H%M%S')
    endtime = st[0].stats.endtime.strftime('%Y%m%dT%H%M%S')
    array_name = f'DISP_{starttime}_{endtime}' # for array name, the time precison is seconds
    if array_name in grp:
        print(f'[WARN] {array_name} alreay exists in {grp}, overwrite!\n')
        h5f.remove_node(grp, array_name)
    ncomp = len(st)
    shape = (ncomp, st[0].stats.npts)
    ca = h5f.create_carray(grp, array_name, atom, shape, filters=filters)
    for i in range(ncomp):
        ca[i, :] = np.array(st[i].data, dtype=np.float32)
    orientation = [(tr.stats.channel[-1], tr.stats.metadata['azimuth'], tr.stats.metadata['dip']) for tr in st]
    dtype = [('code', 'S1'), ('azimuth', float), ('dip', float)]
    ca.attrs['orientation'] = np.array(orientation, dtype=dtype)
    ca.attrs['response_cutoff_frequency'] = [max_flc, min_frc]
    ca.attrs['filter_type'] = 'cosine_sac_taper'
    ca.attrs['filter_limit'] = rmresp_prefilt
    ca.attrs['starttime'] = st[0].stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
    ca.attrs['endtime'] = st[0].stats.endtime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
    ca.attrs['sampling_rate'] = st[0].stats.sampling_rate
    ca.attrs['npts'] = st[0].stats.npts
    ca.attrs['data_type'] = 'DISP'
    ca.attrs['data_unit'] = 'meter'

# save contents in config.yaml
h5f.root._v_attrs['config'] = config

# save event info
h5f.root._v_attrs['event'] = event_data

h5f.flush()
h5f.close()

print('==================================\n')
print(f'[INFO] END: {datetime.datetime.now()}')
