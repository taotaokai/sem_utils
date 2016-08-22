subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_scale_dmodel"
  print '(a)', "    - scale model perturbation "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_scale_dmodel "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <model_name> "
  print '(a)', "    <dmodel_dir> <dmodel_name> "
  print '(a)', "    <max_ratio> "
  print '(a)', "    <out_dir> <out_dmodel_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  max(abs(dmodel_scaled/model)) = max_ratio"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds old model files"
  print '(a)', "  (string) model_name:  proc000***_reg1_<model_name>.bin"
  print '(a)', "  (string) dmodel_dir:  directory holds dmodel files"
  print '(a)', "  (string) dmodel_name:  proc000***_reg1_<dmodel_name>.bin"
  print '(a)', "  (float) max_ratio:  maximum ratio "
  print '(a)', "  (string) out_dir:  out directory for new model"
  print '(a)', "  (string) out_dmodel_name:  proc000***_reg1_<out_dmodel_name>.bin"
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
  integer, parameter :: nargs = 9
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_name
  character(len=MAX_STRING_LEN) :: dmodel_dir
  character(len=MAX_STRING_LEN) :: dmodel_name
  real(dp) :: max_ratio
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_dmodel_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: model
  ! dmodel 
  real(dp), dimension(:,:,:,:), allocatable :: dmodel
  ! scale factor
  real(dp) :: max_ratio_local, max_ratio_all
  real(dp) :: scale_factor

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_add_dmodel_lamda_mu_rho_to_tiso: check your inputs."
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
  read(args(4),'(a)') model_name
  read(args(5),'(a)') dmodel_dir 
  read(args(6),'(a)') dmodel_name
  read(args(7),*) max_ratio
  read(args(8),'(a)') out_dir
  read(args(9),'(a)') out_dmodel_name

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(model(NGLLX,NGLLY,NGLLZ,nspec), &
           dmodel(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== get scale factor
  max_ratio_local = 0.0
  do iproc = myrank, (nproc-1), nrank
    ! read old models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    ! read dmodel
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, dmodel_name, dmodel)
    ! get max ratio
    max_ratio_local = max(max_ratio_local, maxval(abs(dmodel/model)))
  enddo
  call max_all_dp(max_ratio_local, max_ratio_all)

  if (myrank == 0) then
    scale_factor = max_ratio/max_ratio_all
  endif
  call bcast_all_singledp(scale_factor)

  print *, "scale_factor = ", scale_factor

  !====== scale dmodel
  do iproc = myrank, (nproc-1), nrank
    ! read dmodel
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, dmodel_name, dmodel)
    ! write out scaled dmodel 
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_dmodel_name, scale_factor*dmodel)
  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
