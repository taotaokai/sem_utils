subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_gll_value_range"
  print '(a)', "    - print value range of gll model "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_scale_dmodel "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <model_tag> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds model files"
  print '(a)', "  (string) model_tag:  proc000***_reg1_<model_tag>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
end subroutine


program xsem_add_dmodel

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 4
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_tag

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! gll 
  real(dp), dimension(:,:,:,:), allocatable :: gll

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] check your inputs."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') model_dir 
  read(args(4),'(a)') model_tag

  ! get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! intialize arrays 
  allocate(gll(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== get scale factor
  print *, "#iproc min max maxabs"
  do iproc = 0, (nproc-1)
    ! read gll 
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, trim(model_tag), gll)
    ! get max amplitude 
    !print *, iproc, minval(gll), maxval(gll), maxval(abs(gll))
    write(*, "(I5,3E15.5)") iproc, minval(gll), maxval(gll), maxval(abs(gll))
  enddo

end program
