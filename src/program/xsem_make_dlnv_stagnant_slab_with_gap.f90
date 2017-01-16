subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_dlnv_stagnant_slab_with_gap "
  print '(a)', "    - add stagnant slab with gap into SEM model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_dlnv_stagnant_slab_with_gap \"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <origin_lat> <origin_lon> <origin_depth> "
  print '(a)', "    <slab_dip> <slab_strike> <slab_half_width>"
  print '(a)', "    <flat_slab_len> <flat_slab_thickness> <flat_slab_gap_angle> "
  print '(a)', "    <slab_dlnv> <out_tag> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  Slab is a planar slice through the sphere, and " 
  print '(a)', "  MOW is an inverted triangular region inside the slab, and"
  print '(a)', "  forms a body of revolution " 
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (float) origin_lat/lon/depth (deg/km): origin point at the middle of the joint line between the bottom planes of the dipping and flat slabs"
  print '(a)', "  (float) slab_dip/strike:  subducting slab dip,strike angles(deg) and thickness(km)"
  print '(a)', "  (float) slab_half_width:  subducting slab half width along the strike direction(km)"
  print '(a)', "  (float) flat_slab_len(km):  flat slab length (km)"
  print '(a)', "  (float) flat_slab_thickness(km):  flat slab thickness (km)"
  print '(a)', "  (float) slab_gap_angle (deg):  angle of gap in the flat slab, the gap opens from the origin point and increases away from the dipping slab "
  print '(a)', "  (float) slab_dlnv:  velocity perturbation in slab"
  print '(a)', "  (string) out_tag:  output tag (proc*_reg1_<tag>.bin) "
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_add_model_stagnant_slab_with_gap
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 14
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  real(dp) :: origin_lat, origin_lon, origin_depth
  real(dp) :: slab_dip, slab_strike, slab_half_width
  real(dp) :: flat_slab_len, flat_slab_thickness, flat_slab_gap_angle, slab_dlnv
  character(len=MAX_STRING_LEN) :: out_tag, out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  type(sem_mesh_data) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec, nspec
  real(dp), allocatable :: dlnv_gll(:,:,:,:)
  real(dp) :: dlnv

  real(dp) :: xyz(3)
  real(dp), dimension(3) :: origin_vec
  real(dp), dimension(3) :: unit_vec_north, unit_vec_east
  real(dp), dimension(3) :: unit_vec_radial, unit_vec_strike, unit_vec_convergent
  real(dp), dimension(3) :: normal_vec_dipping_slab
  real(dp) :: dipping_slab_thickness
  real(dp) :: dist_along_radial_direction
  real(dp) :: dist_along_strike_direction
  real(dp) :: dist_along_convergent_direction
  real(dp) :: dist_along_dipping_slab_normal
  real(dp) :: xyz_proj_to_flat_slab(3), angle_from_convergent_direction

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
  read(args(6), *) slab_dip
  read(args(7), *) slab_strike
  read(args(8), *) slab_half_width
  read(args(9), *) flat_slab_len
  read(args(10), *) flat_slab_thickness
  read(args(11), *) flat_slab_gap_angle
  read(args(12), *) slab_dlnv
  read(args(13), '(a)') out_tag
  read(args(14), '(a)') out_dir

  !===== geometric parameters of slab
  ! use ECEF coordinate frame

  ! change unit: degree -> radians
  origin_lat = origin_lat * DEGREES_TO_RADIANS
  origin_lon = origin_lon * DEGREES_TO_RADIANS

  slab_dip = slab_dip * DEGREES_TO_RADIANS
  slab_strike = slab_strike * DEGREES_TO_RADIANS
  flat_slab_gap_angle = flat_slab_gap_angle * DEGREES_TO_RADIANS

  ! slab origin point
  call geographic_lla2ecef(origin_lat, origin_lon, -1000.0*origin_depth, &
                           origin_vec(1), origin_vec(2), origin_vec(3))
  origin_vec = origin_vec / R_EARTH
  print *, origin_vec

  ! non-dimenionalize length by R_EARTH_KM
  slab_half_width = slab_half_width/R_EARTH_KM
  flat_slab_len = flat_slab_len/R_EARTH_KM
  flat_slab_thickness = flat_slab_thickness/R_EARTH_KM
  origin_depth = origin_depth/R_EARTH_KM

  dipping_slab_thickness = cos(slab_dip)*flat_slab_thickness

  ! radial, north and east unit vectors at origin point
  unit_vec_radial = origin_vec/(sum(origin_vec**2))**0.5
  unit_vec_north = [ -sin(origin_lat)*cos(origin_lon), -sin(origin_lat)*sin(origin_lon), cos(origin_lat) ]
  unit_vec_east = [ -sin(origin_lon), cos(origin_lon), 0.0_dp ]

  ! unit vector of slab strike
  unit_vec_strike = cos(slab_strike)*unit_vec_north + sin(slab_strike)*unit_vec_east
  ! plate convergence direction at origin point
  unit_vec_convergent = -sin(slab_strike)*unit_vec_north + cos(slab_strike)*unit_vec_east 
  ! normal vector of the dipping slab plane
  normal_vec_dipping_slab = cos(slab_dip)*unit_vec_radial + sin(slab_dip)*unit_vec_convergent

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
            xyz = mesh_data%xyz_glob(:,iglob) - origin_vec

            dist_along_radial_direction = sum(xyz*unit_vec_radial)
            dist_along_strike_direction = sum(xyz*unit_vec_strike)

            dlnv = 0.0

            ! within the slab width and above the flat slab bottom plane
            if (abs(dist_along_strike_direction) <= slab_half_width .and. &
                dist_along_radial_direction >= 0.0 ) then

              dist_along_convergent_direction = sum(xyz * unit_vec_convergent)
              dist_along_dipping_slab_normal = sum(xyz * normal_vec_dipping_slab)

              ! inside dipping slab
              if (dist_along_convergent_direction <= 0 .and. &
                  dist_along_dipping_slab_normal >= 0 .and. &
                  dist_along_dipping_slab_normal <= dipping_slab_thickness) &
              then
                !print *, ispec, "inside dipping slab"
                dlnv = slab_dlnv
              endif

              ! inside flat slab
              if (dist_along_convergent_direction >= 0 .and. &
                  dist_along_convergent_direction <= flat_slab_len .and. &
                  dist_along_radial_direction <= flat_slab_thickness) &
              then
                !print *, ispec, "inside flat slab"
                dlnv = slab_dlnv

                ! project xyz to the flat slab plane
                xyz_proj_to_flat_slab = xyz - sum(xyz*unit_vec_radial)*unit_vec_radial
                ! get angle from convergent direction
                angle_from_convergent_direction = atan2( &
                  sum(xyz_proj_to_flat_slab*unit_vec_strike), &
                  sum(xyz_proj_to_flat_slab*unit_vec_convergent))

                ! inside slab gap
                if (abs(angle_from_convergent_direction) < flat_slab_gap_angle) then
                  !print *, ispec, "inside flat slab gap"
                  dlnv = 0.0
                endif
              endif

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
