subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_vsv2_kappa2_thomsen_to_gll"
  print '(a)', "    - convert model (vsv2,kappa2,eps,gamma,delta) to gll model (vph,vpv,vsv,vsh,eta)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_model_vsv2_kappa2_eps_gamma_to_gll "
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
  print '(a)', "  (string) model_suffix:  read in proc*_reg1_[vsv2,kappa2,eps,gamma,delta]<model_suffix>.bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[vph,vpv,vsv,vsh,eta].bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_model_vsv2_kappa2_eps_gamma_to_gll

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
  real(dp), dimension(:,:,:,:), allocatable :: vph2,vpv2,vsh2,vsv2,eta
  real(dp), dimension(:,:,:,:), allocatable :: kappa2,eps,gamma,delta

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

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  allocate(vph2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(kappa2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(delta(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion,   'vsv2'//trim(model_suffix), vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'kappa2'//trim(model_suffix), kappa2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion,    'eps'//trim(model_suffix), eps)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion,  'gamma'//trim(model_suffix), gamma)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion,  'delta'//trim(model_suffix), delta)

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        eps(:,:,:,ispec) = 0.0
        gamma(:,:,:,ispec) = 0.0
        delta(:,:,:,ispec) = 0.0
      endif
    enddo

    vpv2 = kappa2*vsv2
    vph2 = (1.0 + 2.0*eps) * vpv2 
    vsh2 = (1.0 + 2.0*gamma) * vsv2 
    eta = ((1.0+delta)*kappa2 - 2.0)/((1.0+2.0*eps)*kappa2 - 2.0)

    ! write models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', sqrt(vph2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', sqrt(vpv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', sqrt(vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', sqrt(vsh2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', eta)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
