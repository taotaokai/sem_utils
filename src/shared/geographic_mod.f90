module geographic
! utilies for geographic coordinates conversions 

  implicit none
  private

  ! working precision
  integer, parameter :: wp = kind(0.d0)

  ! constants 
  real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)
  real(wp), parameter :: DEG2RAD = PI / 180.0_wp

  ! WGS84 defining parameters (wikipedia:geodetic_datum)
  real(wp), parameter :: wgs84_a = 6378137.0_wp
  real(wp), parameter :: wgs84_invf = 298.257223563_wp
  ! derived constants 
  real(wp), parameter :: wgs84_f = 1.0 / wgs84_invf
  real(wp), parameter :: wgs84_ecc2 = wgs84_f * (2.0 - wgs84_f)
  real(wp), parameter :: wgs84_one_minus_ecc2 = 1.0 - wgs84_ecc2

  ! exported opertations 
  public :: geodetic2ecef


contains

!///////////////////////////////////////////////////////////////////////////////
subroutine geodetic2ecef(lat, lon, height, x, y, z)
!-convert from geodetic(lat,lon,height) to ECEF(x,y,z)
! reference datum: WGS84 
!
!-inputs: lat, lon (degree), height (meter)
!-outputs: x, y, z(meter)
  
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

end module
