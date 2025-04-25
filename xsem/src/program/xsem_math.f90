subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_math "
  print '(a)', "    -  math operation between two SEM gll models"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_math \"
  print '(a)', "    <nproc> <mesh_dir>"
  print '(a)', "    <model_dir_1> <model_tags_1> "
  print '(a)', "    <model_dir_2> <model_tags_2> "
  print '(a)', "    <math_op>"
  print '(a)', "    <out_dir> <model_tags> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir_1/2:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags_1/2:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (string) math_op:  math operations, e.g. add, sub, mul, div, rdiff, radd"
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. rdiff: relative difference betwee V1 and V2 as (V1-V2)/V2 = V1/V2 - 1"
  print '(a)', "     radd: add relative difference into a reference model (V1 + 1)*V2"
  print '(a)', "  3. divlog: log(V1/V2)"
  print '(a)', "     expmul: exp(V1)*V2"

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
  integer, parameter :: nargs = 9
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir_1, model_tags_1
  character(len=MAX_STRING_LEN) :: model_dir_2, model_tags_2
  character(len=MAX_STRING_LEN) :: math_op
  character(len=MAX_STRING_LEN) :: out_dir, out_model_tags

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank

  ! model names
  integer :: nmodel_1, nmodel_2, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names_1(:)
  character(len=MAX_STRING_LEN), allocatable :: model_names_2(:)
  character(len=MAX_STRING_LEN), allocatable :: out_model_names(:)

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), allocatable :: gll_model_1(:,:,:,:,:), gll_model_2(:,:,:,:,:)

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
  read(args(3), '(a)') model_dir_1
  read(args(4), '(a)') model_tags_1
  read(args(5), '(a)') model_dir_2
  read(args(6), '(a)') model_tags_2
  read(args(7), '(a)') math_op
  read(args(8), '(a)') out_dir
  read(args(9), '(a)') out_model_tags

  !===== parse model tags

  call sem_utils_delimit_string(model_tags_1, ',', model_names_1, nmodel_1)
  call sem_utils_delimit_string(model_tags_2, ',', model_names_2, nmodel_2)
  call sem_utils_delimit_string(out_model_tags, ',', out_model_names, nmodel)

  if (nmodel /= nmodel_1 .or. nmodel /= nmodel_2) then
    if (myrank == 0) then
      print *, '[ERROR] nmodel should be all the same!'
      call abort_mpi() 
    endif
  endif
  call synchronize_all()

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names_1=', (trim(model_names_1(i))//"  ", i=1,nmodel)
    print *, '# model_names_2=', (trim(model_names_2(i))//"  ", i=1,nmodel)
    print *, '# out_model_names=', (trim(out_model_names(i))//"  ", i=1,nmodel)
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
  allocate(gll_model_1(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll_model_2(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    call sem_io_read_gll_file_n(model_dir_1, iproc, iregion, &
                                model_names_1, nmodel, gll_model_1)

    call sem_io_read_gll_file_n(model_dir_2, iproc, iregion, &
                                model_names_2, nmodel, gll_model_2)

    ! math operations
    select case (trim(math_op))
      ! binary operation
      case ('add')
        gll_model_1 = gll_model_1 + gll_model_2
      case ('sub')
        gll_model_1 = gll_model_1 - gll_model_2
      case ('mul')
        gll_model_1 = gll_model_1 * gll_model_2
      case ('div')
        gll_model_1 = gll_model_1 / gll_model_2
      case ('rdiff')
        gll_model_1 = gll_model_1/gll_model_2 - 1.0_dp
      case ('radd')
        gll_model_1 = (gll_model_1 + 1.0_dp) * gll_model_2
      case ('divlog')
        gll_model_1 = log(gll_model_1/gll_model_2)
      case ('expmul')
        gll_model_1 = exp(gll_model_1)*gll_model_2
      case default
        print *, "[ERROR] unrecognized operation: ", trim(math_op)
        stop
    endselect

    ! write out result
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, &
                                out_model_names, nmodel, gll_model_1)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
