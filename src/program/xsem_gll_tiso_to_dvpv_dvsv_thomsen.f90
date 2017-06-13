subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_gll_tiso_to_dvpv_dvsv_eps_gamma_delta"
  print '(a)', "    - convert GLL(vph,vpv,vsv,vsh,eta) model into (dvpv,dvsv,eps,gamma,delta)"
  print '(a)', "      based on 1D reference model (vpv_ref,vsv_ref)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_gll_tiso_to_dvpv_dvsv_eps_gamma \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <out_dir> <flag_backward>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int)    nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vsv_ref,vpv_ref,vph,vpv,vsv,vsh,eta].bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[dvpv,dvsv,eps,gamma,delta].bin"
  print '(a)', "  (int)    flag_backward:  flag forward (0) or backward (1) conversion"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. Thomsen's parameters: esp = (C11/C33 - 1)/2; gamma = (C66/C44 -1)/2; delta = (C13 + 2*C44)/C33 - 1"
  print '(a)', "     delta is different from the origional definition by Thomsen, but same to the first order of delta"
end subroutine


program xsem_gll_tiso_to_dvpv_dvsv_eps_gamma_delta

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
  integer :: flag_backward

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vsv_ref,vpv_ref,vph,vpv,vsh,vsv,eta
  real(dp), dimension(:,:,:,:), allocatable :: dvpv,dvsv,eps,gamma,delta

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
  read(args(5),*) flag_backward

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
  allocate(vpv_ref(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vph(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dvpv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dvsv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(delta(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv_ref', vpv_ref)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv_ref', vsv_ref)

    ! read models
    if (flag_backward == 1) then
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'dvpv', dvpv)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'dvsv', dvsv)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eps', eps)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'gamma', gamma)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'delta', delta)

      ! enforce isotropy for element with ispec_is_tiso = .false.
      do ispec = 1, nspec
        if (.not. mesh_data%ispec_is_tiso(ispec)) then
          eps(:,:,:,ispec) = 0.0
          gamma(:,:,:,ispec) = 0.0
          delta(:,:,:,ispec) = 0.0
        endif
      enddo

      vpv = (1.0 + dvpv)*vpv_ref
      vsv = (1.0 + dvsv)*vsv_ref
      vph = vpv*(1.0 + 2.0*eps)**0.5
      vsh = vsv*(1.0 + 2.0*gamma)**0.5
      eta = ((1.0+delta)*vpv**2 - 2.0*vsv**2) / (vph**2 - 2.0*vsv**2)

      ! write models
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', vpv)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', vsv)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', vph)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', vsh)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', eta)

    else

      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh)
      call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)

      dvpv = vpv/vpv_ref - 1.0
      dvsv = vsv/vsv_ref - 1.0
      eps = 0.5*(vph**2/vpv**2 - 1.0)
      gamma = 0.5*(vsh**2/vsv**2 - 1.0)
      delta = (eta*(vph**2 - 2.0*vsv**2) + 2.0*vsv**2)/vpv**2 - 1.0

      ! enforce isotropy for element with ispec_is_tiso = .false.
      do ispec = 1, nspec
        if (.not. mesh_data%ispec_is_tiso(ispec)) then
          eps(:,:,:,ispec) = 0.0
          gamma(:,:,:,ispec) = 0.0
          delta(:,:,:,ispec) = 0.0
        endif
      enddo

      ! write models
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dvpv', dvpv)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dvsv', dvsv)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eps', eps)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'gamma', gamma)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'delta', delta)
    endif

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
