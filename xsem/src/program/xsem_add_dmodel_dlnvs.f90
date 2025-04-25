subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_dlnvs"
  print '(a)', "    - add model perturbation "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_dlnvs_kappa_thomsen_elliptic"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <dmodel_dir> <dmodel_suffix>"
  print '(a)', "    <step_length> "
  print '(a)', "    <min_dlnvs> <max_dlnvs>"
  print '(a)', "    <output_dmodel> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  new_model = model + step_length * dmodel"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds old model files proc000***_reg1_dlnvs.bin"
  print '(a)', "  (string) dmodel_dir:  directory holds dmodel files"
  print '(a)', "  (string) dmodel_suffix:  proc000***_reg1_dlnvs<dmodel_suffix>.bin"
  print '(a)', "  (float) step_length:  step size multiplied on dmodel"
  print '(a)', "  (float) min_dlnvs:  minimum value"
  print '(a)', "  (float) max_dlnvs:  maximum value"
  print '(a)', "  (int) output_dmodel:  flag out dmodel (0:no, 1:yes) proc000***_reg1_dlnvs_dmodel.bin"
  print '(a)', "  (string) out_dir:  out directory for new model proc000***_reg1_dlnvs.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_add_dmodel

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 10
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: dmodel_dir
  character(len=MAX_STRING_LEN) :: dmodel_suffix
  real(dp) :: step_length
  real(dp) :: min_dlnvs
  real(dp) :: max_dlnvs
  integer :: output_dmodel
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs_new
  ! dmodel 
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs_dmodel

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] wrong argument number."
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
  read(args(4),'(a)') dmodel_dir
  read(args(5),'(a)') dmodel_suffix
  read(args(6),*) step_length
  read(args(7),*) min_dlnvs
  read(args(8),*) max_dlnvs
  read(args(9),*) output_dmodel
  read(args(10),'(a)') out_dir

  ! process input
  if (dmodel_suffix == 'x') then
    dmodel_suffix = ''
  endif
  if ( min_dlnvs>max_dlnvs ) then
    print *, "[ERROR] wrong min/max value range."
    call abort_mpi()
  endif

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(dlnvs(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dlnvs_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dlnvs_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== create new model
  do iproc = myrank, (nproc-1), nrank

    print *, "====== proc ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read old models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'dlnvs', dlnvs)

    ! read dmodel
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'dlnvs'//trim(dmodel_suffix), dlnvs_dmodel)

    ! add dmodel with step_length
    dlnvs_new = dlnvs + step_length*dlnvs_dmodel

    ! limit model value range
    where(dlnvs_new < min_dlnvs) dlnvs_new = min_dlnvs
    where(dlnvs_new > max_dlnvs) dlnvs_new = max_dlnvs

    print *, "dlnvs min/max = ", minval(dlnvs_new), maxval(dlnvs_new)

    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dlnvs', dlnvs_new)

    ! write out dmodel
    if (output_dmodel == 1) then
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dlnvs_dmodel', dlnvs_new - dlnvs)
    endif

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
