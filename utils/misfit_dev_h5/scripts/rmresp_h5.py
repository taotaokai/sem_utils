import sys
import glob
import numpy as np
import scipy
from obspy.signal.invsim import cosine_sac_taper
from obspy import read, read_inventory, Inventory
from obspy.signal.util import _npts2nfft
import tables as pt
from rmresp_utils import get_cutoff_freq
import matplotlib.pyplot as plt

""" data structure
/<net>_<sta>_<loc>_<chan>/<type>_<starttime>_<endtime>[:]
                          attrs: starttime, sampling_rate, npts,
                                 type(e.g. RAW, DISP, VEL),
                                 unit(e.g. count, meter, m/s)
                                 filter

e.g.

/BE_MEM_00_HHE/DISP_20200321T004951_20200321T011950[:]
               attrs: starttime='2020-03-21T00:49:51.008393',
                      sampling_rate=100.0,
                      npts=18000,
                      type='DISP',
                      unit='meter',
                      filter_limit=[0.01, 0.015, 10, 12],
                      filter_type='cosine_sac_taper'

NOTE precision of time is microsecond

"""
h5f = pt.open_file("data_remove_resp.h5", mode="w", title="data")
# taper_alpha = 0.05 # frequency domain taper width [(1-taper_alpha)*flc, flc, frc, (1+taper_alpha)*frc]
sampling_rate = 10 # sampling rate used
lc_taper, hc_taper = 0.01, 0.05 # frequency domian cosine window [lc, lc + lc_taper, hc - hc_taper, hc]
                                # will introduce side lope of duration ~ 1/taper seconds on each side of the spike response

xml_list = glob.glob('stations/*.xml')
mseed_list = glob.glob('waveforms/*.mseed')
# xml_list = glob.glob('test/*.xml')
# mseed_list = glob.glob('test/*.mseed')

# read in all staitonxml files
inv = Inventory()
for xml_file in xml_list:
    try:
        inv1 = read_inventory(xml_file)
    except Exception as err:
        print(f'cannot read {xml_file}, skip!\n{err}')
        continue
    inv.extend(inv1)

# process each mseed file
for mseed_file in mseed_list:
    print(f'[INFO] =========================================')
    print(f'[INFO] processing {mseed_file}')

    try:
        st = read(mseed_file)
    except Exception as err:
        print(f"[WARN] cannot read {mseed_file}, skip.\n{err}")
        continue
    st.merge(-1) # merge continuous segments if possible
    st.detrend('linear')

    for tr in st:
        print(f'[INFO] processing trace {tr}')

        # tr1 = tr.copy()

        # check if response data exist
        try:
            resp = inv.get_response(tr.id, tr.stats.starttime)
        except Exception as err:
            print(f"[Error] cannot get channel response data for {tr.id}, skip.\n{err}")
            continue

        # # raw data: /nslc/raw[npts]
        # shape = (tr.stats.npts,)
        # atom = pt.Atom.from_dtype(tr.data.dtype)
        # filters = pt.Filters(complevel=5, complib='zlib')
        # ca = h5f.create_carray(nslc_g, 'raw', atom, shape, filters=filters)
        # ca[:] = tr.data

        tr.detrend('linear')

        # apply lowpass filter before resampling
        fs = tr.stats.sampling_rate
        npts = tr.stats.npts
        npad = int(2.0 / hc_taper * fs)
        print(f'[INFO] resample npad = {npad}')
        nfft = scipy.fft.next_fast_len(npts + npad)
        freqs = np.fft.rfftfreq(nfft, d=1/fs)

        sig_spectrum = np.fft.rfft(tr.data, nfft)
        nyq_freq = sampling_rate * 0.5
        print(f'[INFO] nyquist frequency = {nyq_freq}')
        resample_hstop = nyq_freq
        resample_hpass = resample_hstop - hc_taper
        resample_prefilt = [-2, -1, resample_hpass, resample_hstop]
        print(f'[INFO] resample_prefilt = {resample_prefilt}')
        taper = cosine_sac_taper(freqs, resample_prefilt) # apply taper near the end
        sig_spectrum *= taper
        tr.data = np.fft.irfft(sig_spectrum, nfft)[:npts]

        # resample by lanczos interpolation
        tr.interpolate(sampling_rate, method='lanczos', a=20)

        # tr.filter('lowpass', freq=resample_hstop, zerophase=True)
        # plt.plot(tr.times(), sig_lp, 'r', tr.times(), tr.data, 'b')
        # plt.show()

        # nfft = tr.stats.npts
        # f = np.fft.rfftfreq(nfft, d=1/sampling_rate)
        # fx = np.fft.rfft(tr.data, nfft)
        # plt.loglog(freqs, np.abs(sig_spectrum), 'k', f, np.abs(fx), 'r')
        # plt.show()

        # get instrument response
        fs = tr.stats.sampling_rate
        # npts = tr.stats.npts
        # nfft1 = _npts2nfft(npts)
        # nfft = scipy.fft.next_fast_len(npts + int(2/flc))
        # print(npts, nfft1, nfft)
        # freqs = np.fft.rfftfreq(nfft, d=1/fs)
        vel_freqs = np.hstack((1.0/np.arange(1000, 0, -1), np.arange(2, fs, 1)))
        try:
            vel_resp = resp.get_evalresp_response_for_frequencies(vel_freqs, output='DEF')
            # disp_resp = resp.get_evalresp_response_for_frequencies(freqs, output='DISP')
        except Exception as err:
            print(f"[Error] failed to evaluate channel response for {tr.id}, skip.\n{err}")
            continue

        # get corner frequency from velocity response spectrum
        try:
            ref_freq = resp.response_stages[0].normalization_frequency
        except Exception as err:
            print(f"[WARN] cannot get response normalization frequency for {tr.id}, set to 1.0 Hz.\n{err}")
            ref_freq = 1.0
        print(f'[INFO] ref_freq = {ref_freq}')

        flc, frc = get_cutoff_freq(vel_freqs, np.abs(vel_resp), ref_freq=ref_freq, cutoff_dB=-3)
        print(f'[INFO] response corner frequency = [{flc}, {frc}]')
        # rmresp_lstop = (1-taper_alpha)*flc
        # rmresp_lpass = flc
        # rmresp_lstop = 0.5 * flc # 0.5 is chosen so that the duration of the filter is ~ 2/flc
        rmresp_lstop = flc
        rmresp_lpass = flc + lc_taper # filter half duration about 1.0/lc_taper sec
        if flc <= 0: # no cutoff on low frequency
            rmresp_lstop = -2
            rmresp_lpass = -1
        # if frc > resample_hstop:
        #     rmresp_hpass, rmresp_hstop = fs, 2*fs # no additional filter needed on high frequency
        #     total_hpass, total_hstop = resample_hpass, resample_hstop
        # else:
        min_hstop = min(resample_hpass, frc)
        rmresp_hstop = min_hstop
        rmresp_hpass = min_hstop - hc_taper
        rmresp_prefilt = [rmresp_lstop, rmresp_lpass, rmresp_hpass, rmresp_hstop]
        print(f'[INFO] rmresp_prefilt = {rmresp_prefilt}')
        if rmresp_lpass >= rmresp_hpass:
            print(f'[WARN] rmresp_lpass({rmresp_lpass}) >= rmresp_hpass({rmresp_hpass}), skip')
            continue

        # apply taper in frequency domain
        npts = tr.stats.npts
        # npad = 100
        # if flc > 0: npad = max(npad, int(2*fs/flc))
        npad = int(2.0 / min(lc_taper, hc_taper) * fs)
        print(f'[INFO] rmresp npad = {npad}')
        nfft = scipy.fft.next_fast_len(npts + npad)
        # nfft = _npts2nfft(npts)
        freqs = np.fft.rfftfreq(nfft, d=1/fs)

        sig_spectrum = np.fft.rfft(tr.data, nfft)
        taper = cosine_sac_taper(freqs, rmresp_prefilt)
        sig_spectrum *= taper

        # remove instrument response
        try:
            disp_resp = resp.get_evalresp_response_for_frequencies(freqs, output='DISP')
        except Exception as err:
            print(f"[Error] failed to evaluate channel response for {tr.id}, skip.\n{err}")
            continue
        inds = (freqs != 0) & (freqs >= rmresp_lstop) & (freqs <= rmresp_hstop)
        sig_spectrum[inds] /= disp_resp[inds]

        # sig_spectrum[1:] /= disp_resp[1:]
        # sig_spectrum[0] = 0.0
        # sig_spectrum[-1] = abs(sig_spectrum[-1]) + 0.0j
        # plt.loglog(freqs, np.abs(sig_spectrum), 'k', freqs, np.abs(fx), 'r')
        # plt.show()

        sig_rmresp = np.fft.irfft(sig_spectrum, nfft)[:npts]

        # create group if not exists: e.g. /BE_MEM_00_HHE/
        nslc = [tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel]
        grp_name = '_'.join(nslc)
        if grp_name not in h5f.root:
            nslc_g = h5f.create_group(h5f.root, grp_name)
            nslc_g._v_attrs['network'] = nslc[0]
            nslc_g._v_attrs['station'] = nslc[1]
            nslc_g._v_attrs['location'] = nslc[2]
            nslc_g._v_attrs['channel'] = nslc[3]
        else:
            print(f'[INFO] {grp_name} already exists!')
            nslc_g = h5f.get_node(h5f.root, grp_name)

        # store in h5 file, e.g. /BE_MEM_00_HHE/DISP_20200321T004951_20200321T011950
        shape = (tr.stats.npts,)
        atom = pt.Atom.from_dtype(np.dtype(np.float32))
        filters = pt.Filters(complevel=5, complib='zlib')
        starttime = tr.stats.starttime.strftime('%Y%m%dT%H%M%S')
        endtime = tr.stats.endtime.strftime('%Y%m%dT%H%M%S')
        array_name = f'DISP_{starttime}_{endtime}' # for array name, the time precison is seconds
        if array_name in nslc_g:
            print(f'[WARN] {array_name} alreay exists in {nslc_g}, skip!')
            continue
        ca = h5f.create_carray(nslc_g, array_name, atom, shape, filters=filters)
        ca[:] = np.array(sig_rmresp, dtype=np.float32)
        ca.attrs['response_corner_frequency'] = [flc, frc]
        ca.attrs['filter_type'] = 'cosine_sac_taper'
        ca.attrs['filter_limit'] = rmresp_prefilt
        ca.attrs['starttime'] = tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f') # microseconds precision
        ca.attrs['sampling_rate'] = tr.stats.sampling_rate
        ca.attrs['npts'] = tr.stats.npts
        ca.attrs['type'] = 'DISP'
        ca.attrs['unit'] = 'meter'

        # tr1.remove_response(inv, water_level=None, taper=False, zero_mean=False, pre_filt=rmresp_prefilt, output='DISP')
        # tr.remove_response(inv, water_level=None, taper=False, zero_mean=False, pre_filt=rmresp_prefilt, output='DISP')
        # # a = np.max(np.abs(sig_rmresp))
        # # a = 1.0
        # # b = np.max(np.abs(tr1.data))
        # plt.plot(tr.times(), sig_rmresp, 'r', tr1.times(), tr1.data, 'k', tr.times(), tr.data, 'b')
        # # plt.plot(tr.times(), sig_rmresp, 'r', tr.times(), tr.data, 'k')
        # print(np.max(np.abs((sig_rmresp - tr.data)) / np.max(np.abs(sig_rmresp))))
        # plt.show()


h5f.flush()
h5f.close()
