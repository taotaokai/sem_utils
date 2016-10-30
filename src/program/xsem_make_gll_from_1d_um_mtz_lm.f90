subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_gll_from_1d_um_mtz_lm"
  print '(a)', "    - make gll files from 1d model: piecewise model(upper-mantle, mantle transition zoen, lower mantle)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_gll_from_1d_um_mtz_lm \"
  print '(a)', "    <nproc> <mesh_dir> \"
  print '(a)', "    <model_um> <model_mtz> <model_lm> "
  print '(a)', "    <out_dir> <out_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_um/mtz/lm:  a text file consisting lines of model points: radius value "
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', "  (string) out_name:  output proc*_reg1_<out_name>.txt"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. model_file must be sorted from small to large radius."
end subroutine


program xsem_make_1d_gll

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_file_um, model_file_mtz, model_file_lm
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec, iglob, igllx, iglly, igllz, ispec
  ! model file
  character(len=MAX_STRING_LEN), allocatable :: model_lines(:)
  integer :: npts_um, npts_mtz, npts_lm
  real(dp), allocatable :: model_radius_um(:), model_value_um(:)
  real(dp), allocatable :: model_radius_mtz(:), model_value_mtz(:)
  real(dp), allocatable :: model_radius_lm(:), model_value_lm(:)
  ! model gll
  real(dp), dimension(:,:,:,:), allocatable :: model_gll
  real(dp) :: R_410km, R_650km, r, r_center, gamma

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] check your inputs."
    stop
  endif

  do i = 1, nargs
      call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_file_um
  read(args(4), '(a)') model_file_mtz
  read(args(5), '(a)') model_file_lm
  read(args(6), '(a)') out_dir
  read(args(7), '(a)') out_name 

  !====== read in model file
  call sem_utils_read_line(model_file_um, model_lines, npts_um)
  allocate(model_radius_um(npts_um), model_value_um(npts_um))
  do i = 1, npts_um
    read(model_lines(i), *) model_radius_um(i), model_value_um(i)
  enddo
  print *, "min/max radius: ", minval(model_radius_um), maxval(model_radius_um)
  print *, "min/max values: ", minval(model_value_um), maxval(model_value_um)

  call sem_utils_read_line(model_file_mtz, model_lines, npts_mtz)
  allocate(model_radius_mtz(npts_mtz), model_value_mtz(npts_mtz))
  do i = 1, npts_mtz
    read(model_lines(i), *) model_radius_mtz(i), model_value_mtz(i)
  enddo
  print *, "min/max radius: ", minval(model_radius_mtz), maxval(model_radius_mtz)
  print *, "min/max values: ", minval(model_value_mtz), maxval(model_value_mtz)

  call sem_utils_read_line(model_file_lm, model_lines, npts_lm)
  allocate(model_radius_lm(npts_lm), model_value_lm(npts_lm))
  do i = 1, npts_lm
    read(model_lines(i), *) model_radius_lm(i), model_value_lm(i)
  enddo
  print *, "min/max radius: ", minval(model_radius_lm), maxval(model_radius_lm)
  print *, "min/max values: ", minval(model_value_lm), maxval(model_value_lm)

  !====== loop each slice
  R_410km = 1.0 - 410.0/6371.0
  R_650km = 1.0 - 650.0/6371.0

  do iproc = 0, (nproc-1)

    print *, "iproc = ", iproc

    ! read mesh data
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec
    if (.not. allocated(model_gll)) then
      allocate(model_gll(NGLLX,NGLLY,NGLLZ,nspec))
    endif

    model_gll = 0.0

    do ispec = 1, nspec

      iglob = mesh_data%ibool(MIDX,MIDY,MIDZ,ispec)
      r_center = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

      ! upper-mantle
      if (r_center > R_410km) then

        do igllz =  1, NGLLZ
          do iglly =  1, NGLLY
            do igllx =  1, NGLLX

              iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
              r = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

              ! linear interpolation 
              do i = 1, npts_um-1
                if ( r >= model_radius_um(i) .and. r < model_radius_um(i+1) ) then
                  exit
                endif
              enddo
              if (r < model_radius_um(1)) then
                model_gll(igllx,iglly,igllz,ispec) = model_value_um(1)
              else if (r > model_radius_um(npts_um)) then
                model_gll(igllx,iglly,igllz,ispec) = model_value_um(npts_um)
              else
                gamma = (r - model_radius_um(i))/(model_radius_um(i+1) - model_radius_um(i))
                model_gll(igllx,iglly,igllz,ispec) = (1.0 - gamma)*model_value_um(i) + gamma*model_value_um(i+1)
              endif

            enddo ! igllx
          enddo ! iglly
        enddo ! igllz

      ! mantle transition zone
      else if (r_center < R_410km .and. r_center > R_650km) then

        do igllz =  1, NGLLZ
          do iglly =  1, NGLLY
            do igllx =  1, NGLLX

              iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
              r = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

              ! linear interpolation 
              do i = 1, npts_mtz-1
                if ( r >= model_radius_mtz(i) .and. r < model_radius_mtz(i+1) ) then
                  exit
                endif
              enddo
              if (r < model_radius_mtz(1)) then
                model_gll(igllx,iglly,igllz,ispec) = model_value_mtz(1)
              else if (r > model_radius_mtz(npts_mtz)) then
                model_gll(igllx,iglly,igllz,ispec) = model_value_mtz(npts_mtz)
              else
                gamma = (r - model_radius_mtz(i))/(model_radius_mtz(i+1) - model_radius_mtz(i))
                model_gll(igllx,iglly,igllz,ispec) = (1.0 - gamma)*model_value_mtz(i) + gamma*model_value_mtz(i+1)
              endif

            enddo ! igllx
          enddo ! iglly
        enddo ! igllz

      ! lower mantle
      else if (r_center < R_650km) then

        do igllz =  1, NGLLZ
          do iglly =  1, NGLLY
            do igllx =  1, NGLLX

              iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
              r = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

              ! linear interpolation 
              do i = 1, npts_lm-1
                if ( r >= model_radius_lm(i) .and. r < model_radius_lm(i+1) ) then
                  exit
                endif
              enddo
              if (r < model_radius_lm(1)) then
                model_gll(igllx,iglly,igllz,ispec) = model_value_lm(1)
              else if (r > model_radius_lm(npts_lm)) then
                model_gll(igllx,iglly,igllz,ispec) = model_value_lm(npts_lm)
              else
                gamma = (r - model_radius_lm(i))/(model_radius_lm(i+1) - model_radius_lm(i))
                model_gll(igllx,iglly,igllz,ispec) = (1.0 - gamma)*model_value_lm(i) + gamma*model_value_lm(i+1)
              endif

            enddo ! igllx
          enddo ! iglly
        enddo ! igllz

      endif

    enddo ! ispec

    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, model_gll)
    print *, "min/max: ", minval(model_gll), maxval(model_gll)

  enddo

end program
