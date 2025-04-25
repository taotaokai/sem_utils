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
# from obspy.clients.filesystem.tsindex import Client
from obspy.clients.filesystem.sds import Client
from obspy.taup import TauPyModel


""" directory structure

EVENT/
     | data/ CMTSOLUTION, stations.txt
     | obs/ net.sta.loc.cha.sac

 stations.txt:
    net|sta|loc|cha,...|lat|lon|ele|dep|az,...|dip,...|freq_win(lc,lp,hp,hc)|resp_win(lc,hc)

    freq_win: filter frequency band, lower cutoff, lower pass, higher pass, higher cutoff
    resp_win: instrument response lower&higher corner frequencies

Example:

C202006211907A/
              | data/ CMTSOLUTION, stations.txt
              | obs/ BE.MEM.00.HH[ENZ].sac

 stations.txt:
    BE|MEM|00|HHE,HHN,HHZ|lat|lon|ele|dep|90,0,0|0,0,-90|0.01,0.02,10,12|0.01,5

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
out_sac_dir = sys.argv[3]
out_station_file = sys.argv[4]

# event
event_data = read_events(CMT_file)[0]
event_id = event_data.event_descriptions[0].text
event_origin = event_data.preferred_origin()

print('==================================\n')
print(f'[INFO] event: {event_origin}\n')

station_fd = open(out_station_file, 'w')

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
        source_depth_in_km=event_origin['depth']/1000.0,
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
        pad_time = 5 # seconds
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
        print(f'[WARN] response lower corner frequency > {config_max_lower_cutoff_freq} Hz, station ignored. ({station_id})\n')
        continue
    if resp_hc < config_min_higher_cutoff_freq:
        print(f'[WARN] response higher corner frequency < {config_min_higher_cutoff_freq} Hz, station ignored. ({station_id})\n')
        continue
    rmresp_lstop = max(resp_lc, config_min_lower_cutoff_freq)
    rmresp_lpass = rmresp_lstop + config_lower_cutoff_taper_width # filter half duration about 1.0/taper_width sec
    nyq_freq = sampling_rate * 0.5
    rmresp_hstop = min(resp_hc, nyq_freq) # anti-aliasing low-pass filter for down-sampling
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

    # write out sac file
    evlo = event_origin['longitude']
    evla = event_origin['latitude']
    evdp = event_origin['depth'] / 1000.0 # km
    stlo = st[0].stats.metadata['latitude']
    stla = st[0].stats.metadata['longitude']
    stel = st[0].stats.metadata['elevation']
    stdp = st[0].stats.metadata['local_depth']
    for tr in st:
        tr.stats.sac = {}
        tr.stats.sac['lcalda'] = True
        tr.stats.sac['kevnm'] = event_id
        tr.stats.sac['evlo'] = evlo
        tr.stats.sac['evla'] = evla
        tr.stats.sac['evdp'] = evdp
        tr.stats.sac['stlo'] = stlo
        tr.stats.sac['stla'] = stla
        tr.stats.sac['stel'] = stel
        tr.stats.sac['stdp'] = stdp
        tr.stats.sac['cmpaz'] = tr.stats.metadata['azimuth']
        tr.stats.sac['cmpinc'] = tr.stats.metadata['dip'] + 90
        sac_file = os.path.join(out_sac_dir, tr.id)
        tr.write(sac_file, format="SAC")

    freq_win = ','.join([f"{x}" for x in rmresp_prefilt])
    resp_win = f'{resp_lc:.3e},{resp_hc}'
    chan_names = ','.join([tr.stats.channel for tr in st])
    chan_azimuths = ','.join(["%.2f"%(tr.stats.metadata['azimuth']) for tr in st])
    chan_dips = ','.join(["%.2f"%(tr.stats.metadata['dip']) for tr in st])

    station_info = f"{network}|{station}|{location}|{chan_names}|{stla}|{stlo}|{stel}|{stdp}|{chan_azimuths}|{chan_dips}|{freq_win}|{resp_win}\n"

    station_fd.write(station_info)
    # ca.attrs['filter_type'] = 'obspy.signal.invsim.cosine_sac_taper'

station_fd.close()

print('==================================\n')
print(f'[INFO] END: {datetime.datetime.now()}')
