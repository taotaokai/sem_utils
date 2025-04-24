subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_gll_to_dlnvs_kappa_thomsen"
  print '(a)', "    - convert GLL model into (dlnvs,kappa,eps,gamma,delta)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_model_gll_to_dlnvs_kappa_thomsen \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vsv_ref,vph,vpv,vsv,vsh,eta].bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[dlnvs,kappa,eps,gamma,delta].bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_model_gll_to_dlnvs_kappa_thomsen_elliptic

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
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vsv_ref,vph,vpv,vsh,vsv,eta
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs,kappa,eps,gamma,delta

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
  allocate(vsv_ref(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vph(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dlnvs(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappa(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(delta(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv_ref', vsv_ref)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)

    dlnvs = log(vsv/vsv_ref)
    kappa = vpv/vsv
    eps = 0.5*(vph**2/vpv**2 - 1.0)
    gamma = 0.5*(vsh**2/vsv**2 - 1.0)
    delta = (eta*(vph**2-2.0*vsv**2) - (vpv**2-2.0*vsv**2))/vpv**2

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        eps(:,:,:,ispec) = 0.0
        gamma(:,:,:,ispec) = 0.0
        delta(:,:,:,ispec) = 0.0
      endif
    enddo

    ! write models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dlnvs', dlnvs)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'kappa', kappa)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eps', eps)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'gamma', gamma)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'delta', delta)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
