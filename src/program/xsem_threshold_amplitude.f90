subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_threshold_amplitude "
  print '(a)', "    - threshold gll file to a given value range"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_threshold_amplitude \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_name> "
  print '(a)', "    <min_value> <max_value> "
  print '(a)', "    <out_dir> <out_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_name>.bin"
  print '(a)', "  (string) model_name:  e.g. mu_kernel "
  print '(a)', "  (float) min_value: minimum value"
  print '(a)', "  (float) max_value: maximum value "
  print '(a)', "  (string) out_dir:  output directory for threshold model files "
  print '(a)', "  (string) out_name: proc*_reg1_<out_name> "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_threshold_amplitude
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables
  ! command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir, model_name
  real(dp) :: min_value, max_value  
  character(len=MAX_STRING_LEN) :: out_dir, out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), allocatable :: model(:,:,:,:)

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] check your input arguments."
      call abort_mpi()
    endif 
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_name
  read(args(5), *) min_value 
  read(args(6), *) max_value
  read(args(7), '(a)') out_dir
  read(args(8), '(a)') out_name

  ! validate input
  if (myrank == 0) then
    if (min_value >= max_value) then
      print *, "[ERROR] min_value must be less than max_value."
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize arrays
  allocate(model(NGLLX,NGLLY,NGLLZ,nspec))

  !===== loop each mesh/model slice
  do iproc = myrank, (nproc-1), nrank
    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    ! threshold 
    where(model < min_value) model = min_value
    where(model > max_value) model = max_value
    print *, "min/max = ", minval(model), maxval(model)
    ! write model
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, model)
  enddo

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
