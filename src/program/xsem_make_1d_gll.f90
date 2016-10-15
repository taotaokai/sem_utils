subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_1d_gll"
  print '(a)', "    - make 1d gll model from the given 1d profile"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_1d_gll \"
  print '(a)', "    <nproc> <mesh_dir> <model_file> <out_dir> <out_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_file:  a text file consisting lines of model points: radius value "
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', "  (string) out_name:  output proc*_reg1_<out_name>.txt"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. model_file must be sorted, I won't do that for you."
end subroutine


program xsem_make_1d_gll

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_file
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
  integer :: npts
  real(dp), allocatable :: model_radius(:), model_value(:)
  ! model gll
  real(dp), dimension(:,:,:,:), allocatable :: model_gll
  real(dp) :: r, gamma

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
  read(args(3), '(a)') model_file
  read(args(4), '(a)') out_dir
  read(args(5), '(a)') out_name 

  !====== read in model file
  call sem_utils_read_line(model_file, model_lines, npts)
  allocate(model_radius(npts), model_value(npts))
  do i = 1, npts
    read(model_lines(i), *) model_radius(i), model_value(i)
  enddo
  print *, "min/max radius: ", minval(model_radius), maxval(model_radius)
  print *, "min/max values: ", minval(model_value), maxval(model_value)

  !====== loop each slice
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
      do igllz =  1, NGLLZ
        do iglly =  1, NGLLY
          do igllx =  1, NGLLX

            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            r = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

            ! linear interpolation 
            do i = 1, npts-1
              if (r>=model_radius(i) .and. r<model_radius(i+1)) then
                exit
              endif
            enddo
            if (r < model_radius(1)) then
              model_gll(igllx,iglly,igllz,ispec) = model_value(1)
            else if (r > model_radius(npts)) then
              model_gll(igllx,iglly,igllz,ispec) = model_value(1)
            else
              gamma = (r - model_radius(i))/(model_radius(i+1) - model_radius(i))
              model_gll(igllx,iglly,igllz,ispec) = (1.0 - gamma)*model_value(i) + gamma*model_value(i+1)
            endif

          enddo ! igllx
        enddo ! iglly
      enddo ! igllz
    enddo ! ispec

    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, model_gll)
    print *, "min/max: ", minval(model_gll), maxval(model_gll)

  enddo

end program
