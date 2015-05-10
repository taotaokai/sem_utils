module geographic
! utilies for geographic coordinates conversions 

  implicit none
  private

  ! working precision: double
  integer, parameter :: dp = kind(0.d0)

  ! constants 
  real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: DEG2RAD = PI / 180.0_dp
  real(dp), parameter :: RAD2DEG = 180.0_dp / PI

  ! WGS84 defining parameters (wikipedia:geodetic_datum)
  real(dp), parameter :: wgs84_a = 6378137.0_dp
  real(dp), parameter :: wgs84_invf = 298.257223563_dp
  ! derived constants 
  real(dp), parameter :: wgs84_f = 1.0_dp / wgs84_invf
  real(dp), parameter :: wgs84_ecc2 = wgs84_f * (2.0_dp - wgs84_f)
  real(dp), parameter :: wgs84_one_minus_ecc2 = 1.0_dp - wgs84_ecc2
  real(dp), parameter :: wgs84_sq_one_minus_f = (1.0_dp - wgs84_f)**2

  ! exported opertations 
  public :: geographic_geodetic2ecef
  public :: geographic_ecef2ned 
  public :: geographic_ecef2latlon0
  public :: rotation_matrix

contains

!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_geodetic2ecef(lat, lon, height, x, y, z)
!-convert from geodetic coordinate(lat,lon,height) to ECEF coordinate(x,y,z)
! ECEF: earth centered earth fixed
! reference datum: WGS84
!
!-inputs: lat,lon[radians], height[meter]
!-outputs: x,y,z[meter]
  
  real(dp), intent(in) :: lat, lon, height
  real(dp), intent(out) :: x, y, z

  real(dp) :: sinlat, coslat, sinlon, coslon
  real(dp) :: normal

  sinlat = sin(lat)
  coslat = cos(lat)
  sinlon = sin(lon)
  coslon = cos(lon)

  normal = wgs84_a / sqrt(1.0_dp - wgs84_ecc2*sinlat**2)

  x = (normal + height) * coslat * coslon
  y = (normal + height) * coslat * sinlon
  z = (normal*wgs84_one_minus_ecc2 + height) * sinlat
    
end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_ecef2ned(x, y, z, olat, olon, oheight, xnorth, yeast, zdown)
!-convert from ECEF(x,y,z) to locat NED(xnorth,yeast,zdown) relative to the
! origin(olat,olon,oheight)
! reference datum: WGS84
!
!-input: x,y,z[meter], olat,olon[radians], oheigth[meter]
!-output: xnorth,yeast,zdown[meter]
  
  real(dp), intent(in) :: x, y, z, olat, olon, oheight
  real(dp), intent(out) :: xnorth, yeast, zdown

  real(dp) :: ox, oy, oz
  real(dp) :: sinlat, coslat, sinlon, coslon
  real(dp) :: dx, dy, dz

  call geographic_geodetic2ecef(olat, olon, oheight, ox, oy, oz)

  sinlat = sin(olat)
  coslat = cos(olat)
  sinlon = sin(olon)
  coslon = cos(olon)

  dx = x - ox
  dy = y - oy
  dz = z - oz

  xnorth = - sinlat * (coslon * dx + sinlon * dy) + coslat * dz
  yeast  = - sinlon * dx + coslon * dy 
  zdown  = - coslat * (coslon * dx + sinlon * dy) - sinlat * dz
    
end subroutine


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
subroutine geographic_ecef2latlon0(x, y, z, lat, lon)
!- get (lat,lon) where the vector (x,y,z) in ECEF intercepts the WGS84 ellipsoid
!
!-input:
! x,y,z: vector(x,y,z), the magnitude is always ignored
!-output:
! lat,lon: geodetic coordinates where the vector(x,y,z) will intercept the
!          ellipsoid

  real(dp), intent(in) :: x, y, z
  real(dp), intent(out) :: lat, lon

  ! get longitude
  lon = atan2(y, x)

  ! convert from geocentric to geolatitude latitude (zero elevation)
  lat = atan2(z, sqrt(x**2 + y**2)/wgs84_sq_one_minus_f)

end subroutine


!///////////////////////////////////////////////////////////////////////////////
end module
