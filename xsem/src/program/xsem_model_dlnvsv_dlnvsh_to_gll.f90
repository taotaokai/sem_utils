subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_dlnvs_to_gll"
  print '(a)', "    - convert model (dlnvs) to gll model (vs) based on a reference vs model "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_model_dlnvs_kappa_thomsen_to_gll"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_suffix> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory containing model files"
  print '(a)', "  (string) model_suffix:  proc*_reg1_[vs_ref,dlnvs]<model_suffix>.bin, use 'x' if no model_suffix present"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_vs.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. vs = vs_ref * exp(dlnvs)"
end subroutine


program xsem_model_dlnvs_to_gll

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_suffix
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
  real(dp), dimension(:,:,:,:), allocatable :: vs_ref
  real(dp), dimension(:,:,:,:), allocatable :: vs
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs

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
  read(args(4),'(a)') model_suffix
  read(args(5),'(a)') out_dir

  ! process inputs
  if (model_suffix == 'x') then
    model_suffix = ''
  endif

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  allocate(vs_ref(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dlnvs(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank

    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vs_ref', vs_ref)

    ! read model perturbations
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'dlnvs'//trim(model_suffix), dlnvs)

    vs = vs_ref * exp(dlnvs)

    ! write models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vs', vs)

  enddo

  !====== finalize
  call synchronize_all()
  call finalize_mpi()

end program
