import os
import sys
import numpy as np
import scipy

from obspy.core import Stream, Trace, Stats
from obspy import read

import tables as pt

def read_h5(h5file):
    with pt.open_file(h5file, mode="a") as h5f:
        data = {}
        stats = Stats()
        stats.sac = {}
        stats.sac['lcalda'] = True
        for evt_g in h5f.root:
            event = evt_g._v_attrs['event']
            stats.sac['evlo'] = evt_g._v_attrs['location']['longitude']
            stats.sac['evla'] = evt_g._v_attrs['location']['latitude']
            stats.sac['evdp'] = evt_g._v_attrs['location']['depth'] / 1000.0 # km
            event_name = event.event_descriptions[0].text
            traces = []
            for sta_g in evt_g:
                stats.network = sta_g._v_attrs['network']
                stats.station = sta_g._v_attrs['station']
                stats.location = sta_g._v_attrs['location']
                stats.sac['stlo'] = sta_g._v_attrs['longitude']
                stats.sac['stla'] = sta_g._v_attrs['latitude']
                stats.sac['stel'] = sta_g._v_attrs['elevation']
                stats.sac['stdp'] = sta_g._v_attrs['local_depth']
                band_inst = sta_g._v_attrs['band_instrument']
                record = sta_g['DATA_DISP']
                ncomp, npts = record.shape
                for icomp in range(ncomp):
                    cmpinf = record._v_attrs['component'][icomp]
                    cmpnm = cmpinf['code']
                    stats.sac['cmpaz'] = cmpinf['azimuth']
                    stats.sac['cmpinc'] = cmpinf['dip'] + 90
                    stats.channel = band_inst + cmpnm.decode()
                    stats.starttime = record._v_attrs['starttime']
                    stats.sampling_rate = record._v_attrs['sampling_rate']
                    # stats.filter_limit = record._v_attrs['filter_param']
                    # stats.type = record.attrs['type']
                    # stats.unit = record.attrs['unit']
                    stats.npts = npts
                    traces.append(Trace(record[icomp, :], stats))
            data[event_name] = Stream(traces=traces)
    return data

if __name__ == '__main__':
    data = read_h5('data_v2.h5')
    sac_dir = 'sac_obs'
    for evt in data:
        for tr in data[evt]:
            path = os.path.join(sac_dir, tr.id)
            tr.write(path, format="SAC")
    # for nslc in set([tr.id[:-1] for tr in st]):
    #     print(nslc)
    #     st.select(id=f'{nslc}?').plot()
