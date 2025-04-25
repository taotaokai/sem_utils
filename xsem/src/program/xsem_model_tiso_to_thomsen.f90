subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_tiso_to_thomsen"
  print '(a)', "    get thomsen parameters(epsilon,delta,gamma) from a vti model (vpv,vph,vsv,vsh,eta)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_model_tiso_to_thomsen \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  Thomsen anisotropy parameters for weak VTI medium (symmetric axis is 3):"
  print '(a)', "    epsilon = (C11 - C33)/C33/2"
  print '(a)', "    gamma   = (C66 - C44)/C44/2"
  print '(a)', "    lamda   = (C13 + C44)/(C33 - C44)"
  print '(a)', "    delta   = (lamda^2 - 1)*(C33 - C44)/C33/2"
  print '(a)', ""
  print '(a)', "  , where C11/rho = Vph^2, C33/rho = Vpv^2, C44/rho = Vsv^2, C66/rho = Vsh^2"
  print '(a)', "          C13 = eta*(C11 - 2*C44)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vpv,vph,vsv,vsh,eta].bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[eps,gamma,delta,eps_minus_delta]"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. for more information: https://en.wikipedia.org/wiki/Transverse_isotropy"

!\begin{align}
!     V_{qP}(\theta) & \approx V_{P0}(1 + \delta \sin^2 \theta \cos^2 \theta + \epsilon \sin^4 \theta) \\
!     V_{qS}(\theta) & \approx V_{S0}\left[1 + \left(\frac{V_{P0}}{ V_{S0}}\right)^2(\epsilon-\delta) \sin^2 \theta \cos^2 \theta\right] \\
!     V_{S}(\theta)  & \approx V_{S0}(1 + \gamma \sin^2 \theta )
!\end{align}
end subroutine


program xsem_model_tiso_to_thomsen

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 4
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: C11, C33, C44, C66, C13, rho
  real(dp), dimension(:,:,:,:), allocatable :: eps, gamma, delta

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

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  allocate(C11(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(C33(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(C44(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(C66(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(C13(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(eps(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(delta(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', C11)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', C33)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', C44)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', C66)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', C13)
    ! convert to squared value
    C11 = C11**2
    C33 = C33**2
    C44 = C44**2
    C66 = C66**2
    C13 = C13*(C11 - 2.0*C44)

    ! thomsen parameters 
    eps   = (C11 - C33)/C33/2.0
    gamma = (C66 - C44)/C44/2.0
    delta = (C13 + C44)/(C33 - C44)
    delta = (delta**2 - 1)*(C33 - C44)/C33/2.0

    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eps', eps)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'gamma', gamma)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'delta', delta)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eps_minus_delta', eps-delta)

  enddo

  !====== Finalize
  if (myrank == 0) close(IOUT)

  call synchronize_all()
  call finalize_mpi()

end program
