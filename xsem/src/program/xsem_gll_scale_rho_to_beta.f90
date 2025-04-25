subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_gll_scale_rho_to_beta"
  print '(a)', "    - scale GLL model (rho) to shear wave perturbation (beta)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_gll_scale_rho_to_beta \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <out_dir> <scaling_factor>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int)    nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[rho0, beta].bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_rho.bin"
  print '(a)', "  (float) scaling_factor:  scales rho to beta"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. rho = rho0 * (1 + scaling_factor*beta)"
end subroutine


program xsem_gll_scale_rho_to_beta

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
  character(len=MAX_STRING_LEN) :: out_dir
  real(dp) :: scaling_factor

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: rho0,beta,rho

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
  read(args(4),'(a)') out_dir
  read(args(5),*) scaling_factor

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  allocate(rho0(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(beta(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho0', rho0)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'beta', beta)

    ! scale vsv perturbations to density
    rho = (1.0 + scaling_factor*beta) * rho0

    ! write models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho', rho)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
