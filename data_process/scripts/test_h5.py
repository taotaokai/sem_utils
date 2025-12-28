import tables as pt

h5 = pt.open_file('test.h5', mode="r")

event = h5.root._v_attrs['event']

origin = event.preferred_origin()

print(origin)

h5.close()
