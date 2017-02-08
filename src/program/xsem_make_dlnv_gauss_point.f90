subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_dlnv_gauss_point "
  print '(a)', "    - creat a gauss point relative perturbation GLL file"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_dlnv_stagnant_slab_with_gap \"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <origin_lat> <origin_lon> <origin_depth> <origin_dlnv> <one_sigma> "
  print '(a)', "    <flag_ellipticity> <out_tag> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "   dlnv = exp(-r^2/sig^2/2), where r is the distance from the origion point"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (float) origin_lat/lon/depth (deg/deg/km): origin point at the center "
  print '(a)', "  (float) origin_dlnv:  relative perturbation at the center"
  print '(a)', "  (float) one_sigma (km): standard deviation "
  print '(a)', "  (int) flag_ellipticity: (0/1) 0: spherical; other: WGS84 ellipsoid"
  print '(a)', "  (string) out_tag:  output tag (proc*_reg1_<tag>.bin) "
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_make_dlnv_gauss_point
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 10
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  real(dp) :: origin_lat, origin_lon, origin_depth
  real(dp) :: origin_dlnv, one_sigma
  integer :: flag_ellipticity
  character(len=MAX_STRING_LEN) :: out_tag, out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  type(sem_mesh_data) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec, nspec
  real(dp), allocatable :: dlnv_gll(:,:,:,:)
  real(dp) :: dlnv

  real(dp) :: origin_xyz(3)
  real(dp) :: xyz(3), r2, one_sigma2

  !====== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_slab_mow_model: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), *) origin_lat
  read(args(4), *) origin_lon
  read(args(5), *) origin_depth
  read(args(6), *) origin_dlnv
  read(args(7), *) one_sigma
  read(args(8), *) flag_ellipticity
  read(args(9), '(a)') out_tag
  read(args(10), '(a)') out_dir

  !===== geometric parameters of slab
  ! use ECEF coordinate frame

  ! change unit: degree -> radians
  origin_lat = origin_lat * DEGREES_TO_RADIANS
  origin_lon = origin_lon * DEGREES_TO_RADIANS

  ! from geodesic to ECEF coordinates, non-dimensionalized
  if (flag_ellipticity /= 0) then
    call geographic_lla2ecef(origin_lat,origin_lon,-1.0*origin_depth*1000.0, origin_xyz(1),origin_xyz(2),origin_xyz(3)) 
    origin_xyz = origin_xyz/R_EARTH
  else
    origin_xyz(1) = (R_EARTH_KM - origin_depth) * cos(origin_lat) * cos(origin_lon)
    origin_xyz(2) = (R_EARTH_KM - origin_depth) * cos(origin_lat) * sin(origin_lon)
    origin_xyz(3) = (R_EARTH_KM - origin_depth) * sin(origin_lat)
    origin_xyz = origin_xyz/R_EARTH_KM
  endif
  
  one_sigma = one_sigma/R_EARTH_KM
  one_sigma2 = one_sigma**2

  !===== loop each mesh slice
  do iproc = 0, (nproc - 1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! allocate model gll array 
    if (allocated(dlnv_gll)) then
      deallocate(dlnv_gll)
    endif
    allocate(dlnv_gll(NGLLX,NGLLY,NGLLZ,nspec))

    ! add slab model on each gll point 
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            ! distance from origin_xyz
            r2 = sum((xyz - origin_xyz)**2)

            ! avoid exp underflow
            if (r2 > 20.0*one_sigma2) then
              dlnv = 0.0_dp
            else
              dlnv = origin_dlnv*exp(-0.5*r2/one_sigma2)
            endif

            ! get dlnv gll
            dlnv_gll(igllx, iglly, igllz, ispec) = dlnv

          enddo
        enddo
      enddo
    enddo

    ! write out model
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_tag, dlnv_gll)

  enddo ! iproc

end program
