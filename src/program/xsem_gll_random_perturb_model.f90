subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_gll_random_perturb_model"
  print '(a)', "    - randomly perturb GLL model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_gll_random_perturb_dvpv_dvsv_thomsen \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_name_list> <min_value> <max_value> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int)    nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_name>.bin"
  print '(a)', "  (string) model_name_list:  comma separated list of model names e.g. alpha,beta,xi"
  print '(a)', "  (float) min_value:  minimum value of random dmodel"
  print '(a)', "  (float) max_value:  maximum value of random dmodel"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_<model_name>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_gll_random_perturb_model

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
  character(len=MAX_STRING_LEN) :: model_name_list
  real(dp) :: min_value
  real(dp) :: max_value
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)
  ! model
  real(dp), dimension(:,:,:,:,:), allocatable :: rand
  real(dp), dimension(:,:,:,:,:), allocatable :: model

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] check your inputs."
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') model_dir 
  read(args(4),'(a)') model_name_list
  read(args(5),*) min_value
  read(args(6),*) max_value 
  read(args(7),'(a)') out_dir

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  call sem_utils_delimit_string(model_name_list, ',', model_names, nmodel)
  allocate(rand(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(model(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

  !====== randomly perturb model

  call RANDOM_SEED()

  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read model
    call sem_io_read_gll_file_n(model_dir, iproc, iregion, model_names, nmodel, model)

    ! randomly perturb model
    call RANDOM_NUMBER(rand)
    rand = min_value + (max_value-min_value)*rand ! restrict to min/max value
    model = model + rand

    ! write out models
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, model_names, nmodel, model)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
