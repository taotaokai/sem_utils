import sys
import pandas as pd
import numpy as np
import pyvista as pv
import simplekml

def vector2latlon_deg(v):
    lon = np.arctan2(v[1], v[0])
    latc = np.arctan2(v[2], (v[0]**2+v[1]**2)**0.5)
    f = 1/299.8
    factor = 1.0 / (1 - f)**2
    lat = np.arctan(factor * np.tan(latc))
    return np.rad2deg(lat), np.rad2deg(lon)

def geodetic_lat2geocentric_lat(geodetic_lat):
    f = 1/299.8
    factor = (1 - f)**2
    return np.arctan(factor * np.tan(geodetic_lat))



sem_parfile = sys.argv[1]

sem_params = pd.read_csv(
    sem_parfile,
    delimiter=r"\s*=\s*",
    header=None,
    comment="#",
    names=["key", "value"],
    dtype=dict(key=object, value=object),
    index_col=["key"],
    engine='python',
).to_dict()["value"]

mesh_central_lat = sem_params['CENTER_LATITUDE_IN_DEGREES']
mesh_central_lon = sem_params['CENTER_LONGITUDE_IN_DEGREES']
mesh_width_xi =    sem_params['ANGULAR_WIDTH_XI_IN_DEGREES']
mesh_width_eta =   sem_params['ANGULAR_WIDTH_ETA_IN_DEGREES']
mesh_gamma_rot =   sem_params['GAMMA_ROTATION_AZIMUTH']

lat_center = float(mesh_central_lat.lower().replace('d', 'e'))
lon_center = float(mesh_central_lon.lower().replace('d', 'e'))
xi_width = float(mesh_width_xi.lower().replace('d', 'e'))
eta_width = float(mesh_width_eta.lower().replace('d', 'e'))
gamma_rot = float(mesh_gamma_rot.lower().replace('d', 'e'))

lat0 = np.deg2rad(lat_center)
lon0 = np.deg2rad(lon_center)

theta = 0.5*np.pi - geodetic_lat2geocentric_lat(lat0)
print(lat0, theta)
phi = lon0

n_center = np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])
v_easting = np.array([-np.sin(phi), np.cos(phi), 0])
v_northing = np.array([-np.cos(theta) * np.cos(phi),
                       -np.cos(theta) * np.sin(phi),
                       np.sin(theta)])

# rotate (v_eath, v_north)
gamma = np.deg2rad(gamma_rot)
v_xi = np.cos(gamma) * v_easting + np.sin(gamma) * v_northing
v_eta = -np.sin(gamma) * v_easting + np.cos(gamma) * v_northing
# print(v_xi, v_eta)

# boundary points
xi = 0.5 * np.deg2rad(xi_width)
eta = 0.5 * np.deg2rad(eta_width)

l_xi = np.tan(xi)
l_eta = np.tan(eta)

# corner points
v = v_xi * l_xi + v_eta * l_eta + n_center
lat0, lon0 = vector2latlon_deg(v)
v = v_xi * l_xi - v_eta * l_eta + n_center
lat1, lon1 = vector2latlon_deg(v)
v = -1 * v_xi * l_xi - v_eta * l_eta + n_center
lat2, lon2 = vector2latlon_deg(v)
v = -1 * v_xi * l_xi + v_eta * l_eta + n_center
lat3, lon3 = vector2latlon_deg(v)

kml = simplekml.Kml()

ls = kml.newlinestring(name="edge", tessellate=1,
                        coords=[(lon0,lat0), (lon1, lat1), (lon2, lat2), (lon3, lat3), (lon0,lat0)])
ls.style.linestyle.width = 5
ls.style.linestyle.color = simplekml.Color.red

pnt = kml.newpoint(name='center', coords=[(lon_center, lat_center)])
pnt.style.iconstyle.scale = 5
pnt.style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/shapes/target.png"

# v1 = np.cross(n_center, v_xi)
# v2 = np.cross(v1, v_xi)

lat4, lon4 = vector2latlon_deg(v_xi)
print(lat4, lon4)
v = n_center * l_xi + v_eta * l_eta + v_xi
lat0, lon0 = vector2latlon_deg(v)
v = n_center * l_xi - v_eta * l_eta + v_xi
lat1, lon1 = vector2latlon_deg(v)
v = -1 * n_center * l_xi - v_eta * l_eta + v_xi
lat2, lon2 = vector2latlon_deg(v)
v = -1 * n_center * l_xi + v_eta * l_eta + v_xi
lat3, lon3 = vector2latlon_deg(v)

ls = kml.newlinestring(name="edge2", tessellate=1,
                        coords=[(lon0,lat0), (lon1, lat1), (lon2, lat2), (lon3, lat3), (lon0,lat0)])
ls.style.linestyle.width = 5
ls.style.linestyle.color = simplekml.Color.blue

pnt = kml.newpoint(name='center2', coords=[(lon4, lat4)])
pnt.style.iconstyle.scale = 5
pnt.style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/shapes/target.png"

kml.save('edge.kml')

# edge points
npts = 100
lines = []
n0 = 0
points = np.zeros((0, 3))

edge_pts = np.zeros((npts, 3))
edge_pts[:] = v_eta * l_eta + n_center
edge_pts += (v_xi * l_xi) * np.linspace(-1, 1, npts).reshape((npts, 1))
edge_pts /= np.sum(edge_pts**2, axis=1, keepdims=True)**0.5
points = np.vstack((points, edge_pts))
line = np.zeros(npts+1, dtype=int)
line[0] = npts
line[1:] = n0 + np.arange(npts)
lines.extend(line)
n0 += npts

edge_pts = np.zeros((npts, 3))
edge_pts[:] = -1 * v_eta * l_eta + n_center
edge_pts += (v_xi * l_xi) * np.linspace(-1, 1, npts).reshape((npts, 1))
edge_pts /= np.sum(edge_pts**2, axis=1, keepdims=True)**0.5
points = np.vstack((points, edge_pts))
line = np.zeros(npts+1, dtype=int)
line[0] = npts
line[1:] = n0 + np.arange(npts)
lines.extend(line)
n0 += npts

edge_pts = np.zeros((npts, 3))
edge_pts[:] = v_xi * l_xi + n_center
edge_pts += (v_eta * l_eta) * np.linspace(-1, 1, npts).reshape((npts, 1))
edge_pts /= np.sum(edge_pts**2, axis=1, keepdims=True)**0.5
points = np.vstack((points, edge_pts))
line = np.zeros(npts+1, dtype=int)
line[0] = npts
line[1:] = n0 + np.arange(npts)
lines.extend(line)
n0 += npts

edge_pts = np.zeros((npts, 3))
edge_pts[:] = -1 * v_xi * l_xi + n_center
edge_pts += (v_eta * l_eta) * np.linspace(-1, 1, npts).reshape((npts, 1))
edge_pts /= np.sum(edge_pts**2, axis=1, keepdims=True)**0.5
points = np.vstack((points, edge_pts))
line = np.zeros(npts+1, dtype=int)
line[0] = npts
line[1:] = n0 + np.arange(npts)
lines.extend(line)
n0 += npts

mesh = pv.PolyData(points, lines=lines)
mesh.save('edge.vtk')
