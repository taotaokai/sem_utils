subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_dmodel_random"
  print '(a)', "    - make random model perturbation "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_dmodel_random "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <min_value> <max_value>"
  print '(a)', "    <out_dir> <out_name> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (float) min_value:  minimum value of random dmodel"
  print '(a)', "  (float) max_value:  maximum value of random dmodel"
  print '(a)', "  (string) out_dir:  out directory for dmodel"
  print '(a)', "  (string) out_name:  proc000***_reg1_<out_name>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_make_dmodel_random

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  real(dp) :: min_value
  real(dp) :: max_value 
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! dmodel
  real(dp), dimension(:,:,:,:), allocatable :: dmodel

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] wrong number of inputs!"
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),*) min_value
  read(args(4),*) max_value 
  read(args(5),'(a)') out_dir
  read(args(6),'(a)') out_name

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(dmodel(NGLLX,NGLLY,NGLLZ,nspec))
 
  call RANDOM_SEED()

  !====== create new model
  do iproc = myrank, (nproc-1), nrank

    ! make random dmodel
    call RANDOM_NUMBER(dmodel)
    ! restrict to min/max value
    dmodel = min_value + (max_value - min_value)*dmodel

    ! write out dmodel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, dmodel)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
