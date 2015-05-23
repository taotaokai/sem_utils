subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xmap_gcircle - get points along a great circle"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xmap_gcircle <lat0> <lon0> <azimuth> <theta_list>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  The great circle is originated at (lat0, lon0) shooting at "
  print '(a)', "    angle <aziumuth>, latitudes and longitudes of the point  "
  print '(a)', "    specified by theta_list is output to stdout"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (float) lat0,lon0: origin point on the great circle (degrees)"
  print '(a)', "  (float) azimuth: shooting azimuth of the great circle (degrees)"
  print '(a)', "  (string) theta_list: list of theta angles (degrees) "
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  1. This auxiliary program is used to plot the great circle path on a map"
  print '(a)', "      in conjuction with xsem_slice_gcircle"
  print '(a)', "    ~~~~~~~~~~~~~~~~"
  print '(a)', "    #!/bin/bash "
  print '(a)', "    ps=<model_tag>.ps "
  print '(a)', "    R=<theta0>/<theta1>/<radius0>/<radius1>"
  print '(a)', "    gmt grdreformat <out_file>?<model_tag> <model_tag>.grd "
  print '(a)', "    gmt grd2cpt <model_tag>.grd -Cseis -R$R -Z -D > cpt "
  print '(a)', "    gmt grdimage <model_tag>.grd -R$R -JPa6iz -Ccpt \"
  print '(a)', "      -BNSEW -Bxa10f5 -Bya100f50 > $ps "
  print '(a)', "    gs $ps" 
  print '(a)', "    ~~~~~~~~~~~~~~~~"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xmap_gcircle

  use sem_constants
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  !-- command line args
  integer, parameter :: nargs = 4
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: theta_list
  real(dp) :: lat0, lon0, azimuth

  !-- local variables
  integer :: i, ier

  !-- theta_list
  integer :: itheta, ntheta
  character(len=MAX_STRING_LEN), allocatable :: theta_lines(:)
  real(dp), allocatable :: theta(:)

  !-- points along great circle 
  real(dp) :: v0(3), v1(3), vnorth(3), veast(3), v_axis(3), vr(3), rotmat(3,3)
  real(dp), allocatable :: lat(:), lon(:)

  !===== read command line arguments

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (len_trim(args(i)) == 0) then
      call selfdoc()
      stop "ERROR: check your input arguments!"
    endif
  enddo

  read(args(1), *) lat0
  read(args(2), *) lon0
  read(args(3), *) azimuth
  read(args(4), '(a)') theta_list

  !===== read theta_list
  call sem_utils_read_line(theta_list, theta_lines, ntheta)
  allocate(theta(ntheta))
  do itheta = 1, ntheta
    read(theta_lines(itheta), *) theta(itheta)
  enddo

  !===== get lat,lon at theta along great circle 

  ! convert units
  lat0 = lat0 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  azimuth = azimuth * DEGREES_TO_RADIANS
  theta = theta * DEGREES_TO_RADIANS

  ! unit direction vector v0 at origin point of the great circle
  call geographic_geodetic2ecef(lat0, lon0, 0.d0, v0(1), v0(2), v0(3))
  v0 = v0 / sqrt(sum(v0**2))

  ! unit direction vector v1 along the shooting azimuth of the great circle
  vnorth = [ - sin(lat0) * cos(lon0), &
             - sin(lat0) * sin(lon0), &
               cos(lat0) ]
  veast = [ - sin(lon0), cos(lon0), 0.d0 ]

  v1 = cos(azimuth) * vnorth + sin(azimuth) * veast

  ! rotation axis = v0 cross-product v1
  v_axis(1) = v0(2)*v1(3) - v0(3)*v1(2)
  v_axis(2) = v0(3)*v1(1) - v0(1)*v1(3)
  v_axis(3) = v0(1)*v1(2) - v0(2)*v1(1)

  ! get vectors at theta(:) in ECEF frame
  allocate(lat(ntheta), lon(ntheta))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    vr = matmul(rotmat, v0)
    call geographic_ecef2latlon0(vr(1), vr(2), vr(3), lat(itheta), lon(itheta))
  enddo

  !===== output results
  lat = lat * RADIANS_TO_DEGREES
  lon = lon * RADIANS_TO_DEGREES
  theta = theta * RADIANS_TO_DEGREES

  do itheta = 1, ntheta
    print '(F10.5,2X,F10.5,2X,F10.5)', theta(itheta), lat(itheta), lon(itheta)
  enddo

end program
