subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh"
  print '(a)', "    - convert GLL (alpha,beta,phi,xi) to (vph,vpv,vsv,vsh) "
  print '(a)', "      based on reference model (vp0,vs0,rho0)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int)    nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[alpha,beta,phi,xi,vp0,vs0].bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[vph,vpv,vsv,vsh].bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. vp0, vs0: voigt average reference models; rho0: reference density model"
  print '(a)', "  3. vp = (1 + alpha)*vp0, vs = (1 + beta)*vs0"
  print '(a)', "     phi = (vph^2 - vpv^2)/vp^2, xi = (vsh^2 - vsv^2)/vs^2"
  print '(a)', "  4. For weak anisotropy (Panning & Romanowicz, 2016) "
  print '(a)', "     vp^2 = 4/5*vph^2 + 1/5*vpv^2, vs^2 = 1/3*vsh^2 + 2/3*vsv^2"
end subroutine


program xsem_gll_alpha_beta_phi_xi_to_vph_vpv_vsv_vsh

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
  ! input model
  real(dp), dimension(:,:,:,:), allocatable :: vp0,vs0
  real(dp), dimension(:,:,:,:), allocatable :: alpha,beta,phi,xi
  ! output model
  real(dp), dimension(:,:,:,:), allocatable :: vp,vs
  real(dp), dimension(:,:,:,:), allocatable :: vph,vpv,vsh,vsv

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
  allocate(vp0(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs0(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(alpha(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(beta(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(phi(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xi(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vp(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vph(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vp0', vp0)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vs0', vs0)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'alpha', alpha)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'beta', beta)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'phi', phi)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'xi', xi)

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        xi(:,:,:,ispec) = 0.0
        phi(:,:,:,ispec) = 0.0
      endif
    enddo

    vp = (1.0 + alpha) * vp0
    vph = sqrt(1.0 + 1.0_dp/5.0_dp*phi) * vp
    vpv = sqrt(1.0 - 4.0_dp/5.0_dp*phi) * vp

    vs = (1.0 + beta) * vs0
    vsv = sqrt(1.0 - 1.0_dp/3.0_dp*xi) * vs
    vsh = sqrt(1.0 + 2.0_dp/3.0_dp*xi) * vs

    ! write models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', vph)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', vsh)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
