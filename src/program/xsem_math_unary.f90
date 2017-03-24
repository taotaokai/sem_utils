subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_math_unary "
  print '(a)', "    -  unary math operation on a single gll file"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_math \"
  print '(a)', "    <nproc> <mesh_dir>"
  print '(a)', "    <model_dir> <model_names> "
  print '(a)', "    <math_op>"
  print '(a)', "    <out_dir> <out_names> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_names:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (string) math_op:  math operations, e.g. abs, zero"
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', "  (string) out_names:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_vertical_slice
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_names
  character(len=MAX_STRING_LEN) :: math_op
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_names

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank
  ! model names
  integer :: nmodel, nmodel_1
  character(len=MAX_STRING_LEN), allocatable :: model_name_list(:)
  character(len=MAX_STRING_LEN), allocatable :: out_name_list(:)
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), allocatable :: model(:,:,:,:,:)

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_math: check your input arguments."
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
  read(args(4), '(a)') model_names
  read(args(5), '(a)') math_op
  read(args(6), '(a)') out_dir
  read(args(7), '(a)') out_names

  !===== parse model tags

  call sem_utils_delimit_string(model_names, ',', model_name_list, nmodel)
  call sem_utils_delimit_string(out_names, ',', out_name_list, nmodel_1)

  if (nmodel /= nmodel_1) then
    if (myrank == 0) then
      print *, '[ERROR] nmodel should be all the same!'
      call abort_mpi()
    endif
  endif
  call synchronize_all()

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_name_list(i))//"  ", i=1,nmodel)
    print *, '# out_names=', (trim(out_name_list(i))//"  ", i=1,nmodel)
  endif

  !===== loop each mesh/model slice
  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! initialize arrays
  allocate(model(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! math operations
    select case (trim(math_op))
      ! binary operation
      case ('abs')
        call sem_io_read_gll_file_n(model_dir, iproc, iregion, model_name_list, nmodel, model)
        model = abs(model)
      case ('zero')
        ! output a gll with zeros, no need to read input gll files
        model = 0.0
      case default
        print *, "[ERROR] unrecognized operation: ", trim(math_op)
        stop
    endselect

    ! write out result
    print *, "min/max = ", minval(model), maxval(model)
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, out_name_list, nmodel, model)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
