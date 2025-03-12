import sys
import os
import glob
import datetime
import numpy as np
import tables as tb
import scipy
import matplotlib.pyplot as plt
import yaml

from obspy.signal.invsim import cosine_sac_taper
from obspy import read, read_inventory, Inventory, read_events, geodetics
# from obspy.signal.util import _npts2nfft
# from obspy.clients.filesystem.tsindex import Client
from obspy.clients.filesystem.sds import Client
from obspy.taup import TauPyModel


""" data structure
/EVENT_ID/ attrs: origin_time, location, event_data
         /NET_STA/ attrs: net,sta,loc,stla,stlo,stel,stdp
                 /DATA_DISP[0:nchan, 0:npts]  (e.g. DATA_RAW, DATA_VEL, SYN_DISP, SYN_VEL)
                       attrs: channel[0:nchan], (e.g. [(name:HHE, azimuth:90, dip:0), ('HHN',...), ('HHZ',...)])
                              starttime, sampling_rate, npts,
                              unit (e.g. count, meter, m/s)
                              filter_type, filter_band
                              response_cutoff_frequency

e.g.

/C202006211907A/ attrs: origin_time='2020-06-21T19:07:00.000000', ...
               /BE_MEM/ attrs: network='BE', station='MEM', location='00',
                      /DATA_DISP[0:3, 0:npts]
                            attrs: channel=[(HHE,90,0), (HHN,0,0), (HHZ,0,-90)],
                                   starttime='2020-03-21T00:49:51.008393',
                                   sampling_rate=10.0,
                                   npts=18000,
                                   unit='meter',
                                   filter_type='cosine_sac_taper'
                                   filter_band=[0.01, 0.02, 10, 12],
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


#
config_yaml = sys.argv[1]
CMT_file = sys.argv[2]
h5_file = sys.argv[3]

# h5 file handle
h5f = tb.open_file(h5_file, mode="a")

# event

event_data = read_events(CMT_file)[0]

event_id = event_data.event_descriptions[0].text
event_origin = event_data.preferred_origin()

print('==================================\n')
print(f'[INFO] event: {event_origin}\n')

# config
with open(config_yaml, 'r') as file:
    config = yaml.safe_load(file)

sampling_rate = config['sampling_rate']

config_zcomp_pattern = config['zcomp_pattern']

config_output_type = config['remove_response']['output_type']
config_output_unit = config['remove_response']['output_unit']
config_ref_resp_type = config['remove_response']['ref_resp_type']
config_ref_freq = config['remove_response']['ref_freq']
config_cutoff_threshold = config['remove_response']['cutoff_threshold']
config_max_lower_cutoff_freq = config['remove_response']['max_lower_cutoff_freq']
config_min_higher_cutoff_freq = config['remove_response']['min_higher_cutoff_freq']
config_min_lower_cutoff_freq = config['remove_response']['min_lower_cutoff_freq']

# filter of frequency domian cosine window [flc, flc + flc_taper, fhc - fhc_taper, fhc]
# response duration ~ 1/min(flc_taper, fhc_taper) seconds
config_lower_cutoff_taper_width = config['remove_response']['lower_cutoff_taper_width']
config_higher_cutoff_taper_width = config['remove_response']['higher_cutoff_taper_width']

time_before_first_arrival = config['time_window']['before_first_arrival_seconds']
# time_after_first_arrival = config['data']['time_window']['after_first_arrival_seconds']
time_after_origin = config['time_window']['after_origin_time_seconds']
config_taper_width = config['time_window']['taper_width']

# for merging traces
misalignment_threshold = config['misalignment_threshold']

stations_path = config['stations_path']
SDS_root = config['SDS_root']

# create /EVENT_ID if not exists:
event_path = f'/{event_id}'
if event_path not in h5f:
    gevent = h5f.create_group('/', event_id)
else:
    gevent = h5f.get_node(event_path)

gevent._v_attrs['event'] = event_data
gevent._v_attrs['origin_time'] = event_origin.time
gevent._v_attrs['location'] = {'longitude': event_origin.longitude,
                               'latitude': event_origin.latitude,
                               'depth': event_origin.depth / 1000.0}
# save configuration
gevent._v_attrs['config'] = config

# taup model
taup_model = TauPyModel(model="ak135")

# SDS client
client = Client(SDS_root)

t0 = event_origin.time - time_before_first_arrival
t1 = event_origin.time + time_after_origin

# query all Z components
# st = client.get_waveforms("*", "*", "*", "[BH]HZ", t0, t1)
st = client.get_waveforms("*", "*", "*", config_zcomp_pattern, t0, t1)

stations = [(tr.stats.network, tr.stats.station) for tr in st]
stations = sorted(set(stations)) # remove duplicates

nslc_list = [(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel) for tr in st]
nslc_list = sorted(set(nslc_list)) # remove duplicates

# in case multiple location/band, keep the one with largest coverage
nslc_filtered = []
for net, sta in stations:
    st1 = st.select(net, sta)
    nslc_list = [(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel) for tr in st1]
    nslc_list = sorted(set(nslc_list)) # remove duplicates
    for i, nslc in enumerate(nslc_list):
        duration = sum([tr.stats.endtime - tr.stats.starttime for tr in st1.select(*nslc)])
        coverage = duration / (t1 - t0)
        nslc_list[i] = (nslc, coverage)
    nslc_list = sorted(nslc_list, key=lambda x:x[1])
    nslc_filtered.append(nslc_list[-1][0])

# process
for network, station, location, channel in nslc_filtered:
    print('==================================')
    station_id = '.'.join([network, station, location, channel[:-1]+'?'])
    print(f'[INFO] {station_id}\n')

    # read StationXML file
    sxml_file = os.path.join(stations_path, f'{network}.{station}.xml')
    try:
        inv = read_inventory(sxml_file)
    except Exception as err:
        print(f"{err}\n[WARN] fail to read {sxml_file}, station ignored. ({station_id})\n")
        continue

    # get channel coordinates
    try:
        seed_id = '.'.join([network, station, location, channel])
        station_coords = inv.get_coordinates(seed_id, event_origin.time)
    except Exception as err:
        print(f"{err}\n[WARN] cannot find station metadata for {seed_id}, station ignored. ({station_id})\n")
        continue

    # calculate event-station distance
    dist, az, baz = geodetics.gps2dist_azimuth(
        event_origin.latitude, event_origin.longitude,
        station_coords['latitude'], station_coords['longitude'])
    dist_degree = geodetics.kilometer2degrees(dist/1000.0)

    # calculate first arrival time using ak135 model
    arrivals = taup_model.get_travel_times(
        source_depth_in_km=event_origin.depth/1000.0,
        distance_in_degree=dist_degree, phase_list=['ttp'])
    first_arrival_time = event_origin.time + min([arr.time for arr in arrivals])

    # data time window
    time_window_starttime = first_arrival_time - time_before_first_arrival
    time_window_endtime = event_origin.time + time_after_origin + config_taper_width
    # # round time to integer seconds
    # if time_window_starttime.microsecond != 0:
    #     time_window_starttime = time_window_starttime.replace(microsecond=0)
    # if time_window_endtime.microsecond != 0:
    #     time_window_endtime = (time_window_endtime + 1).replace(microsecond=0)

    # get waveforms
    try:
        pad_time = 1 # seconds
        st = client.get_waveforms(network, station, location, channel[0:2]+'?',
                                  time_window_starttime - pad_time, time_window_endtime + pad_time) #, merge=False)
    except Exception as err:
        print(f"{err}\n[WARN] failed to get waveforms, station ignored. ({station_id})\n")
        continue

    # merge adjacent traces with misalignment less than sub-sampling interval
    # NOTE merge is done in client.get_waveforms(), which defaults merge=-1
    #      default misalignment_threshold=0.01 should be fine
    # st.merge(-1, misalignment_threshold=misalignment_threshold)

    # ignore traces with gaps or overlaps
    orientation_codes = set([tr.stats.channel[-1] for tr in st]) # get unique orientation codes
    for orientation in orientation_codes:
        st1 = st.select(component=orientation)
        if len(st1) != 1:
            print(f'{st1.__str__(extended=True)}')
            print(f'[WARN] gaps/overlaps found, {st1[0].id} ignored. ({station_id})\n')
            for tr in st1: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace left, station ignored. ({station_id})\n')
        continue

    # check time coverage
    traces_to_remove = []
    for tr in st:
        dt = tr.stats.delta
        if tr.stats.starttime > time_window_starttime:
            print(f'{tr}\n[WARN] trace starttime {tr.stats.starttime} is later than required {time_window_starttime}, {tr.id} ignored. ({station_id})\n')
            traces_to_remove.append(tr)
        elif tr.stats.endtime < time_window_endtime:
            print(f'{tr}\n[WARN] trace endtime {tr.stats.endtime} is earlier than required {time_window_endtime}, {tr.id} ignored. ({station_id})\n')
            traces_to_remove.append(tr)
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace left, station ignored. ({station_id})\n')
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
            print(f"{err}\n[WARN] cannot get instrument response, {tr.id} ignored. ({station_id})\n")
            traces_to_remove.append(tr)
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace left, station ignored. ({station_id})\n')
        continue

    # check if all components have the same coordinates
    for key in ['latitude','longitude','elevation','local_depth']:
        val0 = st[0].stats.metadata[key]
        if not all([tr.stats.metadata[key] == val0 for tr in st[1:]]):
            print(f'[WARN] coordinates are not the same across all components, only one is kept. ({station_id})\n')
            break

    # get response corner frequency on response spectrum
    traces_to_remove = []
    # frequency samples [1/1000, 1/999, ... , 0.5, 1, 1.5, ...] Hz
    ref_freqs = np.hstack((
        1.0/np.arange(1000, 1, -1),
        np.arange(1, sampling_rate/2 +
                  config_higher_cutoff_taper_width+1.5, 0.5)))
    for tr in st:
        fs = tr.stats.sampling_rate
        try:
            # get instrument response spectrum
            resp = tr.stats.response
            ref_resp = resp.get_evalresp_response_for_frequencies(
                    ref_freqs, output=config_ref_resp_type)
        except Exception as err:
            print(f"{err}\n[WARN] failed to evaluate instrument response, {tr.id} ignored. ({station_id})\n")
            traces_to_remove.append(tr)
            continue
        flc, fhc = _get_cutoff_freq(ref_freqs, np.abs(ref_resp), ref_freq=config_ref_freq, cutoff_threshold=config_cutoff_threshold)
        tr.stats.resposne_corner_frequency = [flc, fhc]
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace left, station ignored. ({station_id})\n')
        continue


    # determine common response frequency band for all components
    resp_lc = max([tr.stats.resposne_corner_frequency[0] for tr in st if 'resposne_corner_frequency' in tr.stats])
    resp_hc = min([tr.stats.resposne_corner_frequency[1] for tr in st if 'resposne_corner_frequency' in tr.stats])
    print(f'[INFO] response cutoff frequency = [{resp_lc}, {resp_hc}] for {station_id}\n')
    if resp_lc > config_max_lower_cutoff_freq:
        print(f'[WARN] response lower cutoff frequency > {config_max_lower_cutoff_freq} Hz, station ignored. ({station_id})\n')
        continue
    if resp_hc < config_min_higher_cutoff_freq:
        print(f'[WARN] response higher cutoff frequency < {config_min_higher_cutoff_freq} Hz, station ignored. ({station_id})\n')
        continue
    rmresp_lstop = max(resp_lc, config_min_lower_cutoff_freq)
    rmresp_lpass = rmresp_lstop + config_lower_cutoff_taper_width # filter half duration about 1.0/taper_width sec
    nyq_freq = sampling_rate * 0.5
    rmresp_hstop = min(resp_hc, nyq_freq)
    rmresp_hpass = rmresp_hstop - config_higher_cutoff_taper_width
    if rmresp_lpass >= rmresp_hpass:
        print(f'[WARN] rmresp_lpass({rmresp_lpass}) >= rmresp_rpass({rmresp_hpass}), station ignored. ({station_id})\n')
        continue
    rmresp_prefilt = [rmresp_lstop, rmresp_lpass, rmresp_hpass, rmresp_hstop]
    print(f'[INFO] rmresp_prefilt = {rmresp_prefilt}\n')


    # remove trend and apply taper
    st.detrend('linear')
    st.taper(0.5, max_length=config_taper_width)

    # remove instrument response and resample
    traces_to_remove = []
    resample_starttime = time_window_starttime
    resample_npts = int((time_window_endtime - time_window_starttime)
                        * sampling_rate)
    # nyq_freq = sampling_rate * 0.5
    # resample_hstop = nyq_freq
    # resample_hpass = resample_hstop - config_higher_cutoff_taper_width
    # resample_prefilt = [-2, -1, resample_hpass, resample_hstop]
    # print(f'[INFO] resample_prefilt = {resample_prefilt}\n')
    for tr in st:
        # apply lowpass filter before resampling
        fs = tr.stats.sampling_rate
        npts = tr.stats.npts
        npad = int(2.0 / min(config_lower_cutoff_taper_width, config_higher_cutoff_taper_width) * fs)
        # print(f'[INFO] resample npad = {npad}')
        nfft = scipy.fft.next_fast_len(npts + npad)
        freqs = np.fft.rfftfreq(nfft, d=1/fs)
        sig_spectrum = np.fft.rfft(tr.data, nfft)
        # print(f'[INFO] nyquist frequency = {nyq_freq}')
        # print(f'[INFO] resample_prefilt = {resample_prefilt}')
        taper = cosine_sac_taper(freqs, rmresp_prefilt)
        sig_spectrum *= taper
        # remove instrument response
        try:
            resp = tr.stats.response
            output_resp = resp.get_evalresp_response_for_frequencies(freqs, output=config_output_type)
        except Exception as err:
            print(f"{err}\n[WARN] failed to evaluate instrument response, {tr.id} ignored. ({station_id})\n")
            traces_to_remove.append(tr)
            continue
        inds = (freqs != 0) & (freqs >= rmresp_lstop) & (freqs <= rmresp_hstop)
        sig_spectrum[inds] /= output_resp[inds]
        tr.data = np.fft.irfft(sig_spectrum, nfft)[:npts]
        # resample by lanczos interpolation
        tr.interpolate(sampling_rate,
                       starttime=resample_starttime,
                       npts=resample_npts, method='lanczos', a=20)
    for tr in traces_to_remove: st.remove(tr)
    if len(st) == 0:
        print(f'[WARN] no trace left, station ignored. ({station_id})\n')
        continue

    # st.taper(0.5, max_length=config_taper_width)

    # for debug
    if _DEBUG:
        st1 = st.copy()
        for tr in st1:
            tr.stats.network = 'RAW_' + tr.stats.network

    # for debug
    if _DEBUG:
        (st1 + st).plot(equal_scale=False)

    # check if Z component exists
    if len(st.select(component='Z')) != 1:
        print(f'[WARN] no Z component, station ignored. ({station_id})\n')
        continue

    # check if the dip angle of Z component is -90.0
    tr = st.select(component='Z')[0]
    dip = tr.stats.metadata['dip']
    if dip != -90.0: #FIXME what about 90.0?
        print(f'[WARN] Z component\'s dip angle ({dip}) != -90, station ignored. ({station_id})\n')
        continue

    # check horizontal orientation codes: either 1,2 or E,N
    traces_to_remove = []
    horizontal_orientations = {o for o in orientation_codes if o != 'Z'}
    if len(horizontal_orientations) > 0:
        if horizontal_orientations != {'1', '2'} and horizontal_orientations != {'E', 'N'}:
            print(f'[WARN] horizontal components ({horizontal_orientations}) != (1,2) or (N,E), ignored. ({station_id})\n')
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
                    print(f'[WARN] horizontal components\' dip angles ({dip1}, {dip2}) != 0, ignored. ({station_id})\n')
                    traces_to_remove += [tr1, tr2]
                # if abs((az1 - az2)%180 / 90 - 1) > 1e-5:
                if abs((az1 - az2)%360) < 10:
                    print(f'[WARN] horizontal components\' azimuth angles ({az1}, {az2}) are close to parallel, ignored. ({station_id})\n')
                    traces_to_remove += [tr1, tr2]
    for tr in traces_to_remove: st.remove(tr)

    # check if no trace left
    if len(st) == 0:
        print(f'[WARN] no trace left, station ignored. ({station_id})\n')
        continue

    if config['enforce_3comp'] and len(st) != 3:
        print(f'[WARN] not exactly 3 components, station ignored. ({station_id})\n')
        continue

    # store waveforms in array, e.g. /<EVENT_ID>/BE_MEM/DATA_DISP[0:nchan, 0:npts]
    atom = tb.Atom.from_dtype(np.dtype(np.float32))
    filters = tb.Filters(complevel=3, complib='zlib')

    # create station group
    # start = tr.stats.starttime.strftime('%Y%m%dT%H%M%S')
    # end = tr.stats.endtime.strftime('%Y%m%dT%H%M%S')
    station_name = f'NSL_{network}_{station}_{location}'
    if station_name not in gevent:
        gsta = h5f.create_group(gevent, station_name)
    else:
        gsta = h5f.get_node(gevent, station_name)
    gsta._v_attrs['network'] = st[0].stats.network
    gsta._v_attrs['station'] = st[0].stats.station
    gsta._v_attrs['location'] = st[0].stats.location
    # gsta._v_attrs['band_instrument'] = st[0].stats.channel[0:2]
    gsta._v_attrs['latitude'] = float(st[0].stats.metadata['latitude'])
    gsta._v_attrs['longitude'] = float(st[0].stats.metadata['longitude'])
    gsta._v_attrs['elevation'] = float(st[0].stats.metadata['elevation'])
    gsta._v_attrs['local_depth'] = float(st[0].stats.metadata['local_depth'])

    # write data trace
    trace_name = f'DATA_{config_output_type}' # for array name, the time precison is seconds
    if trace_name in gsta:
        print(f'[WARN] {trace_name} alreay in {gsta}, overwrite. ({station_id})\n')
        h5f.remove_node(gsta, trace_name)
    nchan = len(st)
    shape = (nchan, tr.stats.npts)
    ca = h5f.create_carray(gsta, trace_name, atom, shape, filters=filters)
    for i in range(nchan):
        ca[i, :] = np.array(st[i].data, dtype=np.float32)
    # set attributes
    dtype = [('code', 'S3'), ('azimuth', float), ('dip', float)]
    component = [(tr.stats.channel, tr.stats.metadata['azimuth'], tr.stats.metadata['dip']) for tr in st]
    ca.attrs['component'] = np.array(component, dtype=dtype)
    ca.attrs['starttime'] = st[0].stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
    # ca.attrs['endtime'] = st[0].stats.endtime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
    ca.attrs['sampling_rate'] = st[0].stats.sampling_rate
    ca.attrs['npts'] = st[0].stats.npts
    ca.attrs['unit'] = config_output_unit
    ca.attrs['filter_type'] = 'cosine_sac_taper'
    ca.attrs['filter_param'] = rmresp_prefilt
    ca.attrs['response_cutoff_frequency'] = [resp_lc, resp_hc]

h5f.flush()
h5f.close()

print('==================================\n')
print(f'[INFO] END: {datetime.datetime.now()}')
