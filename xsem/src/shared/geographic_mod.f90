module geographic
! utilies for geographic coordinates conversions 

  implicit none
  private

  ! working precision: double
  integer, parameter :: dp = kind(0.d0)

  ! constants 
  real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
! real(dp), parameter :: DEG2RAD = PI / 180.0_dp
! real(dp), parameter :: RAD2DEG = 180.0_dp / PI

  ! Reference Ellipsoid
  ! WGS84 defining parameters (wikipedia:geodetic_datum)
  real(dp), parameter :: wgs84_a = 6378137.0
  real(dp), parameter :: wgs84_invf = 298.257223563
  ! derived constants 
  real(dp), parameter :: wgs84_f = 1.0/wgs84_invf
  real(dp), parameter :: wgs84_one_minus_f = 1.0 - wgs84_f
  real(dp), parameter :: wgs84_b = wgs84_a * wgs84_one_minus_f
  real(dp), parameter :: wgs84_sq_one_minus_f = wgs84_one_minus_f**2
  real(dp), parameter :: wgs84_e2 = 1.0 - wgs84_sq_one_minus_f
  real(dp), parameter :: wgs84_ep2 = 1.0/wgs84_sq_one_minus_f - 1.0

  ! exported opertations 
  public :: geographic_lla2ecef
  public :: geographic_ecef2ned
  public :: geographic_ecef2lla
  public :: geographic_ecef2ll_zeroalt
  public :: geographic_geocentric_lat
  public :: rotation_matrix

contains

!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_lla2ecef(lat, lon, alt, x, y, z)
!-convert geodetic(lat,lon,alt) to ECEF(x,y,z) coordinates
! ECEF: earth centered earth fixed
! Reference ellipsoid: WGS84
! alt: ellipsoidal heigth
!
!-inputs: lat,lon[rad], alt[meter]
!-outputs: x,y,z[meter]
  
  real(dp), intent(in) :: lat, lon, alt
  real(dp), intent(out) :: x, y, z

  real(dp) :: sinlat, coslat, sinlon, coslon
  real(dp) :: N 

  sinlat = sin(lat)
  coslat = cos(lat)
  sinlon = sin(lon)
  coslon = cos(lon)

  N = wgs84_a / sqrt(1.0_dp - wgs84_e2*sinlat**2)

  x = (N + alt) * coslat * coslon
  y = (N + alt) * coslat * sinlon
  z = (N*wgs84_sq_one_minus_f + alt) * sinlat
    
end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_ecef2ned(x, y, z, lat0, lon0, alt0, north, east, down)
!-convert ECEF(x,y,z) to NED(north,east,down) referred to (lat0,lon0,alt0)
! reference ellipsoid: WGS84
!
!-input: x,y,z[meter], lat0,lon0[rad],alt0[meter]
!-output: xnorth,yeast,zdown[meter]
  
  real(dp), intent(in) :: x, y, z, lat0, lon0, alt0
  real(dp), intent(out) :: north, east, down

  real(dp) :: x0, y0, z0
  real(dp) :: sinlat, coslat, sinlon, coslon
  real(dp) :: dx, dy, dz

  call geographic_lla2ecef(lat0, lon0, alt0, x0, y0, z0)

  sinlat = sin(lat0)
  coslat = cos(lat0)
  sinlon = sin(lon0)
  coslon = cos(lon0)

  dx = x - x0
  dy = y - y0
  dz = z - z0

  north = - sinlat * (coslon * dx + sinlon * dy) + coslat * dz
  east  = - sinlon * dx + coslon * dy 
  down  = - coslat * (coslon * dx + sinlon * dy) - sinlat * dz

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_ecef2lla(x, y, z, lat, lon, alt)
!-convert ECEF(x,y,z) to (lat,lon,alt)
! reference ellipsoid: WGS84
!
!-input: x,y,z[meter]
!-output: lat,lon[rad], alt[meter]
  
  real(dp), intent(in) :: x, y, z
  real(dp), intent(out) :: lat, lon, alt

  real(dp) :: p, theta, N

  lon = atan2(y,x)

  p = sqrt(x**2+y**2)
  theta = atan2(z, p*wgs84_one_minus_f)

  lat = atan2(z + wgs84_ep2*wgs84_b*sin(theta)**3, p - wgs84_e2*wgs84_a*cos(theta)**3)
  N = wgs84_a/sqrt(1.0_dp - wgs84_e2*sin(lat)**2)
  alt = p/cos(lat) - N

  !correct for numerical instability in altitude near exact poles:
  !(after this correction, error is about 2 millimeters, which is about
  !the same as the numerical precision of the overall function)
  if (p < 1.0) then
    alt = z - wgs84_b
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_ecef2ll_zeroalt(x, y, z, lat, lon)
!- get (lat,lon) where the vector ECEF(x,y,z) intercepts the WGS84 ellipsoid
!
!-input:
! x,y,z: vector(x,y,z), the magnitude is always ignored
!-output:
! lat,lon[rad]: geodetic coordinates where the vector(x,y,z) will intercept the
!          ellipsoid

  real(dp), intent(in) :: x, y, z
  real(dp), intent(out) :: lat, lon

  ! get longitude
  lon = atan2(y,x)

  ! convert from geocentric to geodetic latitude
  ! tan(phi_geocentric) = (1-f)^2 * tan(phi_geodetic)
  lat = atan2(z, sqrt(x**2+y**2)*wgs84_sq_one_minus_f)

end subroutine


!///////////////////////////////////////////////////////////////////////////////
real(dp) function geographic_geocentric_lat(lat)
!- convert geodesic latitude to geocentric latitude 
!
!-input:
!   lat: geodesic latitude, radians
!-return:
!   geocentric_lat: radians

  real(dp), intent(in) :: lat

  ! convert from geodetic to geocentric latitude
  ! tan(phi_geocentric) = (1-f)^2 * tan(phi_geodetic)
  geographic_geocentric_lat = atan(wgs84_sq_one_minus_f*tan(lat))

end function

!///////////////////////////////////////////////////////////////////////////////
subroutine rotation_matrix(v_axis, theta, R)
!-rotate matrix with axis (v_axis) by angle (theta)
!
!-input:
! v_axis: vector of rotation axis
! theta: rotation angle (radians)
!
!-output:
! R: rotation matrix

  real(dp), intent(in) :: v_axis(3), theta
  real(dp), intent(out) :: R(3,3)

  real(dp) :: sint, cost, one_minus_cost, v(3)

  sint = sin(theta)
  cost = cos(theta)
  one_minus_cost = 1.0 - cost

  ! normalize v_axis to unit vector
  v = v_axis / sqrt(sum(v_axis**2))

  ! rotation matrix
  R(1,1) = cost + one_minus_cost * v_axis(1)**2
  R(2,2) = cost + one_minus_cost * v_axis(2)**2
  R(3,3) = cost + one_minus_cost * v_axis(3)**2

  R(1,2) = one_minus_cost*v_axis(1)*v_axis(2) - sint * v_axis(3)
  R(2,1) = one_minus_cost*v_axis(1)*v_axis(2) + sint * v_axis(3)

  R(1,3) = one_minus_cost*v_axis(1)*v_axis(3) + sint * v_axis(2)
  R(3,1) = one_minus_cost*v_axis(1)*v_axis(3) - sint * v_axis(2)

  R(2,3) = one_minus_cost*v_axis(2)*v_axis(3) - sint * v_axis(1)
  R(3,2) = one_minus_cost*v_axis(2)*v_axis(3) + sint * v_axis(1)

end subroutine


!///////////////////////////////////////////////////////////////////////////////
end module
