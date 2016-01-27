subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_diff_relative "
  print '(a)', "    -  math operation between two SEM gll models"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_math \"
  print '(a)', "    <nproc> <mesh_dir>"
  print '(a)', "    <model_dir_1> <model_tags_1> "
  print '(a)', "    <model_dir_2> <model_tags_2> "
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
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir_1, model_tags_1
  character(len=MAX_STRING_LEN) :: model_dir_2, model_tags_2

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank

  ! model names
  integer :: nmodel, nmodel_2
  character(len=MAX_STRING_LEN), allocatable :: model_names_1(:)
  character(len=MAX_STRING_LEN), allocatable :: model_names_2(:)

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

  !===== parse model tags

  call sem_utils_delimit_string(model_tags_1, ',', model_names_1, nmodel)
  call sem_utils_delimit_string(model_tags_2, ',', model_names_2, nmodel_2)

  if (nmodel /= nmodel_2) then
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

    ! read two model files
    call sem_io_read_gll_file_n(model_dir_1, iproc, iregion, &
                                model_names_1, nmodel, gll_model_1)

    call sem_io_read_gll_file_n(model_dir_2, iproc, iregion, &
                                model_names_2, nmodel, gll_model_2)

    ! relative difference
    gll_model_1 = (gll_model_1 - gll_model_2) / gll_model_1

    ! write out result
    print '("rank= ",I3," min_dlnV= ",E10.3," max_dlnV= ",E10.3)', iproc, &
      minval(gll_model_1), maxval(gll_model_1)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
