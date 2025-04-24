subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_model_geographic_region"
  print '(a)', "    - add model perturbation in a geographic region"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_model_geographic_region \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <lat0> <lat1> <lon0> <lon1> <r0> <r1>"
  print '(a)', "    <dlnV> <flag_output_dlnV> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0/1:  latitude range [lat0,lat1] (degree) "
  print '(a)', "  (float) lon0/1:  longitude range [lon0,lon1] (degree) "
  print '(a)', "  (float) r0/1:  radius range [r0, r1] (km)"
  print '(a)', "  (float) dlnV:  velocity perturbation at plume axis"
  print '(a)', "  (int) flag_output_dlnV:  (1/0) 1=output dlnV, 0=don't output dlnV"
  print '(a)', "  (string) out_dir:  output directory of new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_add_model_horizontal_layer
  
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
  integer :: nproc, flag_output_dlnV
  real(dp) :: lat0, lat1, lon0, lon1, r0, r1, dlnV0

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

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
  real(dp) :: xyz(3), lat, lon, r
  real(dp) :: dlnV
  real(dp), allocatable :: dlnV_gll(:,:,:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_geographic_region: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), *) lat0
  read(args(6), *) lat1
  read(args(7), *) lon0
  read(args(8), *) lon1
  read(args(9), *) r0
  read(args(10), *) r1
  read(args(11), *) dlnV0
  read(args(12), *) flag_output_dlnV
  read(args(13), '(a)') out_dir

  ! validate input arguments
  if (flag_output_dlnV /= 0 .and. flag_output_dlnV /= 1) then
    print *, "[ERROR] flag_output_dlnV must be 0 or 1"
    stop
  endif

  !===== geometric parameters of plume

  ! change unit: degree -> radians
  lat0 = lat0 * DEGREES_TO_RADIANS
  lat1 = lat1 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  lon1 = lon1 * DEGREES_TO_RADIANS

  ! normalize length by EARTH_R_KM
  r0 = r0 / EARTH_R_KM
  r1 = r1 / EARTH_R_KM

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

    if (flag_output_dlnV == 1) then
      if (allocated(dlnV_gll)) then
        deallocate(dlnV_gll)
      endif
      allocate(dlnV_gll(NGLLX,NGLLY,NGLLZ,nspec))
    endif

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
                                model_names, nmodel, model_gll)

    ! add slab model on each gll point
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            ! get (lat, lon, r) from xyz
            lat = atan2(xyz(3), sqrt(xyz(1)**2 + xyz(2)**2))
            lon = atan2(xyz(2), xyz(1))
            r = sqrt(sum(xyz**2))

            ! get the model perturbation
            dlnV = 0.0_dp

            if (lat >= lat0 .and. lat <= lat1 .and. &
                lon >= lon0 .and. lon <= lon1 .and. &
                  r >= r0   .and.   r <= r1 ) then
              dlnV = dlnV0
            endif

            ! get the new model
            model_gll(:,igllx,iglly,igllz,ispec) = (1.0_dp + dlnV) * &
              model_gll(:, igllx,iglly,igllz,ispec)

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
