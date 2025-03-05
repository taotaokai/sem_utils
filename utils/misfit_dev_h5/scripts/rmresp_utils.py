import sys
import numpy as np
import scipy
import obspy
from obspy import Stream, read, read_inventory
from obspy.signal.invsim import cosine_sac_taper


def get_cutoff_freq(freqs, amp, ref_freq=-1, cutoff_dB=-3):
    if ref_freq > 0:
        ref_idx = int(np.interp(ref_freq, freqs, np.arange(len(freqs))))
    else:
        ref_idx = np.argmax(amp)

    ref_amp = amp[ref_idx]
    threshold_amp = ref_amp * 10**(cutoff_dB/20.0)

    # initial cutoff frequency
    left_cutoff_freq, right_cutoff_freq = -10, np.max(freqs)

    inds = np.flatnonzero(amp[:ref_idx] <= threshold_amp)
    if inds.size > 0:
        left_cutoff_idx = np.max(inds)
        left_cutoff_freq = freqs[left_cutoff_idx]

    inds = np.flatnonzero(amp[ref_idx:] <= threshold_amp)
    if inds.size > 0:
        right_cutoff_idx = ref_idx + np.min(inds)
        right_cutoff_freq = freqs[right_cutoff_idx]

    return left_cutoff_freq, right_cutoff_freq


# get corner frequencies of instrument response
def get_resp_cutoff_freq(resp, cutoff_dB=-3, delta_freq=0.001):
    assert(type(resp) is obspy.core.inventory.response.Response)

    sampling_rate = -1.0
    for stage in resp.response_stages[::-1]:
        if (stage.decimation_input_sample_rate is not None
            and stage.decimation_factor is not None):
            sampling_rate = (stage.decimation_input_sample_rate / stage.decimation_factor)
            break
    if sampling_rate < 0.0:
        print('[WARN] unable to determine sampling rate from response, default to 100.0 Hz!')
        sampling_rate = 100.0

    nfreq = int(0.5 * sampling_rate / delta_freq)
    freqs = np.arange(nfreq) * delta_freq
    try:
        cmplx_response = resp.get_evalresp_response_for_frequencies(freqs, output='DEF')
        amp = np.abs(cmplx_response)
    except Exception as err:
        raise(f"[ERROR] cannot eval response for {resp}")

    ref_freq = -1.0
    try:
        ref_freq = resp.response_stages[0].normalization_frequency
    except Exception as err:
        print(f"[WARN] cannot get response normalization frequency from {resp}, set to 1.0 Hz.\n{err}")
        ref_freq = 1.0
    if ref_freq <= 0: ref_freq = 1.0
    print(f"[INFO] response normalization_frequency {ref_freq}")

    return get_cutoff_freq(freqs, amp, ref_freq=ref_freq, cutoff_dB=cutoff_dB)


# remove instrument response
def remove_response(stream, inv, output='DISP', taper_alpha=0.01):
    """
    stream: obspy.core.stream.Stream
    inventory: obspy.core.inventory.Inventory
    """
    assert(type(stream) is obspy.core.stream.Stream)
    assert(type(inv) is obspy.core.inventory.inventory.Inventory)

    output_traces = []
    for tr in stream:
        print('\n=========================================')
        print('[INFO] processing trace: ', tr)

        fs = tr.stats.sampling_rate
        npts = tr.stats.npts
        nfft = scipy.fft.next_fast_len(2 * npts)
        freqs = np.fft.rfftfreq(nfft, d=1/fs)
        try:
            resp = inv.get_response(tr.id, tr.stats.starttime)
            disp_resp = resp.get_evalresp_response_for_frequencies(freqs, output=output)
        except Exception as err:
            print(f"[Error] cannot get channel response data for {tr.id}, skip.\n{err}")
            continue

        # get corner frequency from velocity spectrum
        try:
            vel_resp = resp.get_evalresp_response_for_frequencies(freqs, output='DEF')
            ref_freq = resp.response_stages[0].normalization_frequency
        except Exception as err:
            print(f"[WARN] cannot get response normalization frequency for {tr.id}, set to 1.0 Hz.\n{err}")
            ref_freq = 1.0
        print(f'[INFO] ref_freq = {ref_freq}')
        flc, frc = get_cutoff_freq(freqs, np.abs(vel_resp), ref_freq=ref_freq, cutoff_dB=-3)
        print(f'[INFO] response corner frequency: flc, frc = {flc}, {frc}')
        # resp_flimit = [flc, (1+taper_alpha)*flc, (1-taper_alpha)*frc, frc]
        # if flc <= 0: # no cutoff in lower end
        #     resp_flimit = [-2, -1, (1-taper_alpha)*frc, frc]
        # print(f'[INFO] resp_flimit = {resp_flimit}')

        # apply taper in frequency domain, and remove response
        sig = scipy.signal.detrend(tr.data, type='constant')
        sig_spectrum = np.fft.rfft(sig, nfft)
        # taper = cosine_sac_taper(freqs, resp_flimit)
        # sig_spectrum *= taper
        inds = (freqs < flc) | (freqs > frc)
        sig_spectrum[inds] = 0
        inds = (freqs != 0) & (freqs >= flc) & (freqs <= frc)
        sig_spectrum[inds] /= disp_resp[inds]
        sig_rmresp = np.fft.irfft(sig_spectrum, nfft)[:npts]

        tr = tr.copy()
        tr.data = sig_rmresp
        # tr.stats.filter = {'cosine_sac_taper': resp_flimit}
        tr.stats.filter = [flc, frc]
        tr.stats.type = output
        output_traces.append(tr)

    return Stream(traces=output_traces)


def remove_response_to_disp_obspy(stream, inv, taper_alpha=0.01):
    assert(type(stream) is obspy.core.stream.Stream)
    assert(type(inv) is obspy.core.inventory.inventory.Inventory)

    output_traces = []
    for tr in stream:
        print('\n=========================================')
        print('[INFO] processing trace: ', tr)

        # get corner frequency from velocity spectrum
        try:
            resp = inv.get_response(tr.id, tr.stats.starttime)
            flc, frc = get_resp_cutoff_freq(resp, cutoff_dB=-3, delta_freq=0.0001)
        except Exception as err:
            print(f"[WARN] cannot get response corner frequency for {tr.id}, skip.\n{err}")
            continue
        print(f'[INFO] response corner frequency: flc, frc = {flc}, {frc}')

        flimit = [flc, (1+taper_alpha)*flc, (1-taper_alpha)*frc, frc]
        print(f'[INFO] remove response pre_filt: {flimit}')
        tr1 = tr.copy()
        tr1.detrend('linear')
        tr1.remove_response(inv, water_level=None, taper=False, pre_filt=flimit, zero_mean=False, output='DISP')
        tr1.stats.flimit = flimit
        tr1.stats.type = 'DISP'
        # tr1.stats.unit = 'meter'
        output_traces.append(tr1)
    return Stream(traces=output_traces)

if __name__ == '__main__':
    st = read('data.mseed')
    inv = read_inventory('station.xml')

    st1 = remove_response_to_disp_obspy(st.select(id='NL.HGN.02.HH?'), inv, taper_alpha=0.001)
    print(st1[0].stats)
    st1.plot()

    st2 = remove_response(st.select(id='NL.HGN.02.HH?'), inv)
    print(st2[0].stats)
    st2.plot()
