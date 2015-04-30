module geographic
! utilies for geographic coordinates conversions 

  implicit none
  private

  ! working precision
  integer, parameter :: wp = kind(0.d0)

  ! constants 
  real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)
  real(wp), parameter :: DEG2RAD = PI / 180.0_wp
  real(wp), parameter :: RAD2DEG = 180.0_wp / PI

  ! WGS84 defining parameters (wikipedia:geodetic_datum)
  real(wp), parameter :: wgs84_a = 6378137.0_wp
  real(wp), parameter :: wgs84_invf = 298.257223563_wp
  ! derived constants 
  real(wp), parameter :: wgs84_f = 1.0 / wgs84_invf
  real(wp), parameter :: wgs84_ecc2 = wgs84_f * (2.0 - wgs84_f)
  real(wp), parameter :: wgs84_one_minus_ecc2 = 1.0 - wgs84_ecc2

  ! exported opertations 
  public :: geographic_geodetic2ecef
  public :: geographic_ecef2ned 
  public :: rotation_matrix

contains

!///////////////////////////////////////////////////////////////////////////////
subroutine geographic_geodetic2ecef(lat, lon, height, x, y, z)
!-convert from geodetic coordinate(lat,lon,height) to ECEF coordinate(x,y,z)
! ECEF: earth centered earth fixed
! reference datum: WGS84
!
!-inputs: lat,lon[degree], height[meter]
!-outputs: x,y,z[meter]
  
  real(wp), intent(in) :: lat, lon, height
  real(wp), intent(out) :: x, y, z

  real(wp) :: sinlat, coslat, sinlon, coslon
  real(wp) :: normal

  sinlat = sin(lat * DEG2RAD)
  coslat = cos(lat * DEG2RAD)
  sinlon = sin(lon * DEG2RAD)
  coslon = cos(lon * DEG2RAD)

  normal = wgs84_a / sqrt(1.0 - wgs84_ecc2*sinlat**2)

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
!-input: x,y,z[meter], olat,olon[degree],oheigth[meter]
!-output: xnorth,yeast,zdown[meter]
  
  real(wp), intent(in) :: x, y, z, olat, olon, oheight
  real(wp), intent(out) :: xnorth, yeast, zdown

  real(wp) :: ox, oy, oz
  real(wp) :: sinlat, coslat, sinlon, coslon
  real(wp) :: dx, dy, dz

  call geographic_geodetic2ecef(olat, olon, oheight, ox, oy, oz)

  sinlat = sin(olat * DEG2RAD)
  coslat = cos(olat * DEG2RAD)
  sinlon = sin(olon * DEG2RAD)
  coslon = cos(olon * DEG2RAD)

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
!

  real(wp), intent(in) :: v_axis(3), theta
  real(wp), intent(out) :: R(3,3)

  real(wp) :: sint, cost, one_minus_cost, v(3)

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
