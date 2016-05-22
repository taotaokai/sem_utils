subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_model_llsvp "
  print '(a)', "    - add a dome-like llsvp anomaly above the CMB"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_model_llsvp \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <llsvp_lat> <llsvp_lon> "
  print '(a)', "    <llsvp_dlnV0> <llsvp_sigma> "
  print '(a)', "    <flag_ellipticity> <flag_output_dlnV> "
  print '(a)', "    <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  The llsvp is a vertical cylinder and the velocity "
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
  print '(a)', "  (float) llsvp_lat/lon:  surface location of llsvp axis(degrees) "
  print '(a)', "  (float) llsvp_dlnV0:  velocity perturbation at llsvp center"
  print '(a)', "  (float) llsvp_sigma:  guassian width of llsvp (km)"
  print '(a)', "  (int) flag_ellipticity: (0/1) 0: spherical; 1: WGS84 ellipsoid"
  print '(a)', "  (int) flag_output_dlnV: (0/1) 1=output dlnV, 0=don't output dlnV"
  print '(a)', "  (string) out_dir:  output directory of new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_add_model_llsvp
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 11
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_dir
  integer :: nproc, flag_ellipticity, flag_output_dlnV
  real(dp) :: llsvp_lat, llsvp_lon, llsvp_dlnV0, llsvp_sigma 

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! llsvp geometry
  real(dp) :: llsvp_v(3), llsvp_xyz0(3)
  real(dp), parameter :: r_cmb = 1.0 - 2891.0/R_EARTH_KM

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
  real(dp) :: xyz(3), r, r_llsvp
  real(dp) :: dlnV
  real(dp), allocatable :: dlnV_gll(:,:,:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_llsvp: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), *) llsvp_lat
  read(args(6), *) llsvp_lon
  read(args(7), *) llsvp_dlnV0
  read(args(8), *) llsvp_sigma
  read(args(9), *) flag_ellipticity
  read(args(10), *) flag_output_dlnV
  read(args(11), '(a)') out_dir

  !-- validate input arguments
  if (flag_ellipticity /= 0 .and. flag_ellipticity /= 1) then
    print *, "[ERROR] flag_ellipticity must be 0/1"
    stop
  endif

  if (flag_output_dlnV /= 0 .and. flag_output_dlnV /= 1) then
    print *, "[ERROR] flag_output_dlnV must be 0 or 1"
    stop
  endif

  !===== geometric parameters of llsvp

  ! change unit: degree -> radians
  llsvp_lat = llsvp_lat * DEGREES_TO_RADIANS
  llsvp_lon = llsvp_lon * DEGREES_TO_RADIANS

  ! normalize length by R_EARTH_KM
  llsvp_sigma = llsvp_sigma / R_EARTH_KM

  ! unit directional vector of slab origin
  if (flag_ellipticity == 0) then
    llsvp_v(1) = cos(llsvp_lat) * cos(llsvp_lon)
    llsvp_v(2) = cos(llsvp_lat) * sin(llsvp_lon)
    llsvp_v(3) = sin(llsvp_lat)
  else
    call geographic_lla2ecef(llsvp_lat, llsvp_lon, 0.0_dp, &
      llsvp_v(1), llsvp_v(2), llsvp_v(3)) 
  endif

  llsvp_v = llsvp_v / sqrt(sum(llsvp_v**2))
  llsvp_xyz0 = r_cmb * llsvp_v

  print *, "# llsvp_v=", llsvp_v
  print *, "# llsvp_xyz0=", llsvp_xyz0

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
            dlnV = llsvp_dlnV0 * &
              exp(-0.5 * sum((xyz-llsvp_xyz0)**2)/llsvp_sigma**2)

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
