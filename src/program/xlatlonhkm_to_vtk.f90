subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xlatlonhkm_to_vtk - convert from geodetic to vtk file"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xlatlonh_to_xyz <latlonh_list> <out_file>"
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
  integer, parameter :: nargs = 2
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: latlonh_list
  character(len=MAX_STRING_LEN) :: out_file

  !-- local variables
  integer :: i, ier

  !-- latlonh_list
  integer :: ipoint, npoint
  character(len=MAX_STRING_LEN), allocatable :: latlonh_lines(:)
  real(dp), allocatable :: lat(:), lon(:), height(:)

  !-- ECEF coordinates 
  real(dp), allocatable :: x(:), y(:), z(:)

  !===== read command line arguments

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (len_trim(args(i)) == 0) then
      call selfdoc()
      stop "ERROR: check your input arguments!"
    endif
  enddo

  read(args(1), '(a)') latlonh_list
  read(args(2), '(a)') out_file

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
  do ipoint = 1, npoint
    call geographic_geodetic2ecef(lat(ipoint), lon(ipoint), height(ipoint), &
                                    x(ipoint),   y(ipoint),      z(ipoint))
  enddo

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
  write(IOUT, '(a,2X,a,2X,a)') "Command: xlatlonh_to_vtk ", &
    trim(latlonh_list), trim(out_file)
  write(IOUT, '(a)') "ASCII"
  write(IOUT, '(a)') "DATASET UNSTRUCTURED_GRID"
  write(IOUT, '(a,2X,I10,2X,a)') "POINTS", npoint, "float"
  ! vtk data
  do ipoint = 1, npoint
    write(IOUT,'(2X,E14.6,2X,E14.6,2X,E14.6)') x(ipoint), y(ipoint), z(ipoint)
  enddo

  close(IOUT)

end program
