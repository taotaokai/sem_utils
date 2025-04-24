subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_gll_scale"
  print '(a)', "    - scale gll model by multiplying a constant factor"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_scale_dmodel "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <model_tag> "
  print '(a)', "    <scale_factor> "
  print '(a)', "    <out_dir> <out_tag> "
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
  print '(a)', "  (float) scale_factor:  scale factor"
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', "  (string) out_tag:  <out_dir>/proc000***_reg1_<out_tag>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
end subroutine


program xsem_gll_scale

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
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_tag
  real(dp) :: scale_factor
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_tag

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
  read(args(5),*) scale_factor
  read(args(6),'(a)') out_dir 
  read(args(7),'(a)') out_tag

  ! get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! intialize arrays 
  allocate(gll(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== scale gll model 
  do iproc = 0, (nproc-1)
    ! read in gll 
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, trim(model_tag), gll)
    ! scale gll 
    gll = scale_factor * gll
    ! output gll
    print *, "afer scale: min/max = ", minval(gll), maxval(gll)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, trim(out_tag), gll)
  enddo

end program
