subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_model_plume "
  print '(a)', "    - add vertical plume model into SEM model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_model_plume \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <plume_lat> <plume_lon> "
  print '(a)', "    <plume_r0> <plume_r1> "
  print '(a)', "    <plume_dlnV0> <plume_sigma> "
  print '(a)', "    <flag_ellipticity> <flag_output_dlnV> "
  print '(a)', "    <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  The plume is a vertical cylinder and the velocity "
  print '(a)', "  perturbation decays as a Gaussian function of the "
  print '(a)', "  distance(d) from the vertical axis: "
  print '(a)', "    V'(d) = V(d) * (1 + dlnV0 * exp(- d^2/2/sigma^2)) "
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) plume_lat/lon:  surface location of plume axis(degrees) "
  print '(a)', "  (float) plume_r0/r1:  radius range [r0,r1] of plume (km) "
  print '(a)', "  (float) plume_dlnV0:  velocity perturbation at plume axis"
  print '(a)', "  (float) plume_sigma:  guassian width of plume (km)"
  print '(a)', "  (int) flag_ellipticity: (0/1) 0: spherical; 1: WGS84 ellipsoid"
  print '(a)', "  (int) flag_output_dlnV: (0/1) 1=output dlnV, 0=don't output dlnV"
  print '(a)', "  (string) out_dir:  output directory of new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_vertical_slice
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 13
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_dir
  integer :: nproc, flag_ellipticity, flag_output_dlnV
  real(dp) :: plume_lat, plume_lon, plume_r0, plume_r1, &
              plume_dlnV0, plume_sigma 

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! plume geometry
  real(dp) :: plume_v(3)

  ! model names
  integer :: nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec, nspec
  ! model
  real(dp), allocatable :: model_gll(:,:,:,:,:)

  !-- calculate model perturbation
  real(dp) :: xyz(3), r, r_plume
  real(dp) :: dlnV
  real(dp), allocatable :: dlnV_gll(:,:,:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_plume: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), *) plume_lat
  read(args(6), *) plume_lon
  read(args(7), *) plume_r0
  read(args(8), *) plume_r1
  read(args(9), *) plume_dlnV0
  read(args(10), *) plume_sigma
  read(args(11), *) flag_ellipticity
  read(args(12), *) flag_output_dlnV
  read(args(13), '(a)') out_dir

  !-- validate input arguments
  if (flag_ellipticity /= 0 .and. flag_ellipticity /= 1) then
    print *, "[ERROR] flag_ellipticity must be 0/1"
    stop
  endif

  if (flag_output_dlnV /= 0 .and. flag_output_dlnV /= 1) then
    print *, "[ERROR] flag_output_dlnV must be 0 or 1"
    stop
  endif

  !===== geometric parameters of plume

  ! change unit: degree -> radians
  plume_lat = plume_lat * DEGREES_TO_RADIANS
  plume_lon = plume_lon * DEGREES_TO_RADIANS

  ! normalize length by R_EARTH_KM
  plume_r0 = plume_r0 / R_EARTH_KM
  plume_r1 = plume_r1 / R_EARTH_KM
  plume_sigma = plume_sigma / R_EARTH_KM

  ! unit directional vector of slab origin
  if (flag_ellipticity == 0) then
    plume_v(1) = cos(plume_lat) * cos(plume_lon)
    plume_v(2) = cos(plume_lat) * sin(plume_lon)
    plume_v(3) = sin(plume_lat)
  else
    call geographic_geodetic2ecef(plume_lat, plume_lon, 0.0_dp, &
      plume_v(1), plume_v(2), plume_v(3)) 
  endif

  plume_v = plume_v / sqrt(sum(plume_v**2))

  print *, "# plume_v=", plume_v

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel) 

  !===== loop each mesh slice
  do iproc = 0, (nproc - 1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! read model gll
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
                                model_names, nmodel, model_gll)

    ! initialize dlnV array
    if (flag_output_dlnV == 1) then
      if (allocated(dlnV_gll)) then
        deallocate(dlnV_gll)
      endif
      allocate(dlnV_gll(NGLLX,NGLLY,NGLLZ,nspec))
    endif

    ! add slab model on each gll point
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            ! get the model perturbation
            dlnV = 0.0_dp

            r = sqrt(sum(xyz**2))
            ! projection on plume axis
            r_plume = sum(xyz*plume_v)

            if (r_plume > 0 .and. r >= plume_r0 .and. r <= plume_r1) then
              dlnV = plume_dlnV0 * &
                exp(-0.5 * (r**2 - r_plume**2) / plume_sigma**2)
            endif

            ! get the new model
            model_gll(:,igllx,iglly,igllz,ispec) = &
              (1.0_dp + dlnV) * model_gll(:, igllx,iglly,igllz,ispec)

            if (flag_output_dlnV == 1) then
              dlnV_gll(igllx,iglly,igllz,ispec) = dlnV
            endif

          enddo
        enddo
      enddo
    enddo

    ! write out model
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, &
      model_names, nmodel, model_gll)

    if (flag_output_dlnV == 1) then
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, "dlnV", dlnV_gll)
    endif

  enddo ! iproc

end program
