import numpy as np
from obspy.signal.invsim import cosine_sac_taper
from obspy import read, read_inventory, UTCDateTime
import matplotlib.pyplot as plt

station_list = [
        ('DK.STE05', 'DK.STE05..HHN'),
        # ('NO.BRBA', 'NO.BRBA.00.HHZ'), # CMG-3E
        # ('IM.EKB', 'IM.EKB..HHZ'), # CMG-3T, accelerator
        # ('FR.AJAC', 'FR.AJAC.00.HHZ'),  # CMG-3ESPC
        # ('IV.MDI', 'IV.MDI..HHZ'), # NANOMETRICS TRILLIUM-40S
        # ('PL.KSP', 'PL.KSP..HHZ'), # STS-2
        # ('NO.NC405', 'NO.NC405.00.BHZ'), # CMT-3T
        ]

for ns, seed_id in station_list:

    inv = read_inventory(f'stations/{ns}.xml')
    resp = inv.get_response(seed_id, '2020-03-21T00:47:09.250501')

    # resp.plot(0.0001, output='VEL')

    freqs = np.logspace(-3, 2, 100)

    r = resp.get_evalresp_response_for_frequencies(freqs, output='VEL')
    max_amp = np.max(np.abs(r))
    amp_norm =  np.abs(r) / max_amp

    # threshold = 10**(-3/20)
    threshold = 0.05
    print(threshold)

    # ind0 = np.min(np.flatnonzero(amp_norm > threshold))
    ind0 = np.max(np.flatnonzero(amp_norm < threshold))
    flc = freqs[ind0]
    print(flc)

    # if lc < 0.01: lc = 0.01
    flc_taper = max(0.01, flc)
    flimit = [flc, flc+flc_taper, 20, 30]

    taper = cosine_sac_taper(freqs, flimit)

    # taper = np.ones_like(freqs)
    # ind = freqs <= lc
    # taper[ind] = (freqs[ind] / lc)**4

    plt.loglog(freqs, amp_norm, 'k', freqs, taper, 'r')
    plt.show()

    plt.loglog(freqs, taper/amp_norm)
    plt.show()
