subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_get_radius"
  print '(a)', "    - get radius of gll points"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_get_radius \"
  print '(a)', "    <iproc> <mesh_dir> <model_dir> <model_names>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_names>.bin"
  print '(a)', "  (string) model_names:  vsv,vpv,gamma,eps"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_get_radius

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 4
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: iproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_names

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i
  ! model_names
  character(len=MAX_STRING_LEN), allocatable :: model_name_list(:)
  integer :: nmodel
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec, iglob, igllx, iglly, igllz, ispec
  ! model gll
  real(dp), dimension(:,:,:,:,:), allocatable :: model_gll
  real(dp) :: r, r_center

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] check your inputs."
    stop
  endif

  do i = 1, nargs
      call get_command_argument(i, args(i))
  enddo
  read(args(1), *) iproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_names

  ! read mesh data 
  call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

  call sem_utils_delimit_string(model_names, ',', model_name_list, nmodel)

  nspec = mesh_data%nspec
  allocate(model_gll(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

  call sem_io_read_gll_file_n(model_dir, iproc, iregion, model_name_list, nmodel, model_gll)

  i = 1 ! avoid fake error report from -fbounds-check in gfortran compiler 
  print *, "# radius r_elem_center ", (/ (trim(model_name_list(i))//"  ", i=1,nmodel) /)
  do ispec = 1, nspec
    do igllz =  1, NGLLZ
      do iglly =  1, NGLLY
        do igllx =  1, NGLLX

          iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
          r = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

          iglob = mesh_data%ibool(MIDX,MIDY,MIDZ,ispec)
          r_center = sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))

          write(*,"(*(E18.8))"), r, r_center, model_gll(:,igllx,iglly,igllz,ispec)

        enddo ! igllx
      enddo ! iglly
    enddo ! igllz
  enddo ! ispec

end program
