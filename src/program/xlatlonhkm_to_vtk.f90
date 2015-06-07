subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xlatlonhkm_to_vtk - convert from geodetic to vtk file"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xlatlonhkm_to_xyz <latlonh_list> <flag_ellipticity> <out_file>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  convert from geodetic coordinates (lat,lon,height) to ECEF (x,y,z)"
  print '(a)', "    and stored as an unstructrured VTK ascii file. The value "
  print '(a)', "    is non-dimensionalized by nominal earth radius (6371.0 km)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) latlonh_list: list of (lat, lon, height) in degrees/km"
  print '(a)', "  (int) flag_ellipticity: (0/1) 0: spherical; 1: WGS84 ellipsoid"
  print '(a)', "  (string) out_file: output vtk file name"
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xlatlonhkm_to_vtk

  use sem_constants
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  !-- command line args
  integer, parameter :: nargs = 3
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: latlonh_list
  integer :: flag_ellipticity
  character(len=MAX_STRING_LEN) :: out_file

  !-- local variables
  integer :: i, ier

  !-- latlonh_list
  integer :: ipoint, npoint
  character(len=MAX_STRING_LEN), allocatable :: latlonh_lines(:)
  real(dp), allocatable :: lat(:), lon(:), height(:)

  !-- ECEF coordinates 
  real(dp) :: r, cos_lat_r
  real(dp), allocatable :: x(:), y(:), z(:)

  !===== read command line arguments

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (len_trim(args(i)) == 0) then
      call selfdoc()
      stop "[ERROR] check your input arguments!"
    endif
  enddo

  read(args(1), '(a)') latlonh_list
  read(args(2), *) flag_ellipticity
  read(args(3), '(a)') out_file

  !-- validate input arguments
  if (flag_ellipticity /= 0 .and. flag_ellipticity /= 1) then
    print *, "[ERROR] flag_ellipticity must be 0/1"
    stop
  endif

  !===== read theta_list
  call sem_utils_read_line(latlonh_list, latlonh_lines, npoint)
  allocate(lat(npoint), lon(npoint), height(npoint))
  do ipoint = 1, npoint
    read(latlonh_lines(ipoint), *) lat(ipoint), lon(ipoint), height(ipoint)
  enddo

  !===== convert from geodetic coordinate to ECEF 

  ! convert units: degree -> radians, km -> meter
  lat = lat * DEGREES_TO_RADIANS
  lon = lon * DEGREES_TO_RADIANS
  height = height * 1000.d0

  ! unit direction vector v0 at origin point of the great circle
  allocate(x(npoint), y(npoint), z(npoint))

  if (flag_ellipticity == 1) then
    do ipoint = 1, npoint
      call geographic_geodetic2ecef(lat(ipoint), lon(ipoint), height(ipoint), &
                                      x(ipoint),   y(ipoint),      z(ipoint))
    enddo
  elseif (flag_ellipticity == 0) then
    do ipoint = 1, npoint
      r = R_EARTH + height(ipoint)
      cos_lat_r = cos(lat(ipoint)) * r
      x(ipoint) = cos(lon(ipoint)) * cos_lat_r
      y(ipoint) = sin(lon(ipoint)) * cos_lat_r
      z(ipoint) = sin(lat(ipoint)) * r 
    enddo
  endif

  ! non-dimensinalize ECEF xyz
  x = x / R_EARTH
  y = y / R_EARTH
  z = z / R_EARTH

  !===== output unstructrued VTK ascii file 
  open(IOUT, file=trim(out_file), status='unknown', &
    form='formatted', action='write', iostat=ier)

  if (ier /= 0) then
    write(*,*) '[ERROR] xlatlonh_to_vtk: failed to open file ', trim(out_file)
    stop
  endif

  ! vtk file header
  write(IOUT, '(a)') "# vtk DataFile Version 2.0"
  write(IOUT, '(a,2X,a,2X,I1,2X,a)') "Command: xlatlonhkm_to_vtk ", &
    trim(latlonh_list), flag_ellipticity, trim(out_file)
  write(IOUT, '(a)') "ASCII"
  write(IOUT, '(a)') "DATASET UNSTRUCTURED_GRID"
  write(IOUT, '(a,2X,I10,2X,a)') "POINTS", npoint, "float"
  ! vtk data
  do ipoint = 1, npoint
    write(IOUT,'(2X,E14.6,2X,E14.6,2X,E14.6)') x(ipoint), y(ipoint), z(ipoint)
  enddo

  close(IOUT)

end program
