import numpy as np
import scipy
from obspy import Stream, read, read_inventory
from obspy.signal.invsim import cosine_sac_taper
import matplotlib.pyplot as plt

fs = 10.0
tlen = 1000
lc, hc = 1/360, 4.5
lc_taper, hc_taper = 0.01, 0.1
# flimit = [0.01, 0.02,  10, 20]
npts = int(fs * tlen)
x = scipy.signal.unit_impulse(npts, idx='mid')

nfft = scipy.fft.next_fast_len(2 * npts)
freqs = np.fft.rfftfreq(nfft, d=1/fs)
Fx = np.fft.rfft(x, nfft)

times = np.arange(npts) / fs

# for alpha in [0.99, 0.95, 0.9, 0.85, 0.8]:
#     flimit = [alpha*lc, lc, 10, 20]
#     Ftaper = cosine_sac_taper(freqs, flimit)
#     Fxfilt = Ftaper * Fx
#     xfilt = np.fft.irfft(Fxfilt, nfft)[:npts]
#     plt.plot(times, xfilt,  c=colors[i], label=f'{alpha}')
#     i += 1

# for alpha in [1.01, 1.05, 1.1, 1.15, 1.2]:
#     flimit = [-2, -1, hc, alpha*hc]
#     Ftaper = cosine_sac_taper(freqs, flimit)
#     Fxfilt = Ftaper * Fx
#     xfilt = np.fft.irfft(Fxfilt, nfft)[:npts]
#     plt.plot(times, xfilt,  c=colors[i], label=f'{alpha}')
#     i += 1

# alpha_list = [0.05, 0.1, 0.2, 0.5, 0.8]
# colors = plt.cm.rainbow(np.linspace(0, 1, len(alpha_list)))
# i = 0
# for alpha in alpha_list:
#     # flimit = [(1-alpha)*lc, lc, hc, (1+alpha)*hc]
#     # flimit = [(1-alpha)*lc, lc, 10, 20]
#     flimit = [lc, lc / (1-alpha), 10, 20]
#     # flimit = [-2, -1, (1-alpha)*hc, hc]
#     Ftaper = cosine_sac_taper(freqs, flimit)
#     Fxfilt = Ftaper * Fx
#     xfilt = np.fft.irfft(Fxfilt, nfft)[:npts]
#     plt.plot(times, xfilt,  c=colors[i], label=f'{alpha}')
#     i += 1

# df_list = [0.005, 0.01, 0.02]
# colors = plt.cm.rainbow(np.linspace(0, 1, len(df_list)))
# i = 0
# for df in df_list:
#     # flimit = [(1-alpha)*lc, lc, hc, (1+alpha)*hc]
#     # flimit = [(1-alpha)*lc, lc, 10, 20]
#     # flimit = [lc, lc + df, 10, 20]
#     flimit = [-2, -1, lc, lc + df]
#     Ftaper = cosine_sac_taper(freqs, flimit)
#     Fxfilt = Ftaper * Fx
#     xfilt = np.fft.irfft(Fxfilt, nfft)[:npts]
#     plt.plot(times, xfilt,  c=colors[i], label=f'{df}')
#     i += 1

lc_list = [20, 50, 100, 200]
flc_list = [0, 0.005, 0.01]
# hc_list = [1, 2, 3, 4]
colors = plt.cm.rainbow(np.linspace(0, 1, len(lc_list)))
flc_taper = 0.01
i = 0
for flc in flc_list:
    # flc = 1/lc
    flimit = [flc, flc+flc_taper, 10, 20]

    # taper = np.ones_like(freqs)
    # ind = freqs <= flc
    # taper[ind] = (freqs[ind] / flc)**4
    # Fxfilt = taper * Fx
    # x1 = np.fft.irfft(Fxfilt, nfft)[:npts]

    # taper = np.ones_like(freqs)
    # ind = freqs <= flc
    # taper[ind] = 0
    # Fxfilt = taper * Fx
    # x2 = np.fft.irfft(Fxfilt, nfft)[:npts]

    # flimit = [0, flc, 20, 30]
    taper = cosine_sac_taper(freqs, flimit)
    Fxfilt = taper * Fx
    x3 = np.fft.irfft(Fxfilt, nfft)[:npts]

    plt.plot(times, x3,  c=colors[i], label=f'flc={flc} Hz')
    # plt.plot(times, x2, '--', c=colors[i])
    # plt.plot(times, x3, '.-', c=colors[i])

    i += 1

# Fxfilt = Fx.copy()
# inds = (freqs < lc) | (freqs > hc)
# # inds = (freqs < lc)
# inds = (freqs > hc)
# Fxfilt[inds] = 0
# xfilt = np.fft.irfft(Fxfilt, nfft)[:npts]
# # plt.plot(times, xfilt, c='k', label=f'boxcar')

# plt.plot(times, x/np.max(x), 'k')
plt.legend()
plt.show()
