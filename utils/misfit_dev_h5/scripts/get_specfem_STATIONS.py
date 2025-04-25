import tables as pt

h5 = pt.open_file('data.h5', mode="r")

stations = [(g._v_attrs['network'],
             g._v_attrs['station'],
             g._v_attrs['latitude'],
             g._v_attrs['longitude'],
             g._v_attrs['elevation'],
             g._v_attrs['local_depth'])
            for g in h5.root]

h5.close()

# stations = sorted(stations, key=lambda x: x[0])
stations.sort()

with open('STATIONS', 'w') as f:
    for sta in stations:
        f.write(f'{sta[1]:8s} {sta[0]:5s} {sta[2]:12f} {sta[3]:12f} {sta[4]:12f} {sta[5]:12f}\n')
