import sys
import numpy as np
import scipy

from obspy.core import Stream, Trace, Stats

import tables as pt

def read_h5(h5file):
    stats = Stats()
    traces = []
    with pt.open_file(h5file, mode="r") as h5f:
        for nslc_g in h5f.root:
            try:
                stats.network = nslc_g._v_attrs['network']
                stats.station = nslc_g._v_attrs['station']
                stats.location = nslc_g._v_attrs['location']
                stats.channel = nslc_g._v_attrs['channel']
            except Exception as err:
                print(f"[Error] cannot get channel attributes for {nslc_g}, skip.\n{err}")
                continue
            for record in nslc_g:
                try:
                    stats.starttime = record.attrs['starttime']
                    stats.sampling_rate = record.attrs['sampling_rate']
                    stats.filter_limit = record.attrs['filter_limit']
                    stats.type = record.attrs['type']
                    stats.unit = record.attrs['unit']
                    stats.npts = len(record[:])
                    traces.append(Trace(record[:], stats))
                except Exception as err:
                    print(f"[Error] cannot get trace data for {tr}, skip.\n{err}")
                    continue
    return Stream(traces=traces)

if __name__ == '__main__':
    st = read_h5('data_remove_resp.h5')
    print(st)
    for nslc in set([tr.id[:-1] for tr in st]):
        print(nslc)
        st.select(id=f'{nslc}?').plot()
