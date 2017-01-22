subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_dlnvs_kappa_thomsen_to_gll_serial"
  print '(a)', "    - convert model (dlnvs,kappa,eps,gamma) to gll model (vph,vpv,vsv,vsh,eta)"
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
  print '(a)', "  (string) model_suffix:  proc*_reg1_[vsv_ref,dlnvs,kappa,eps,gamma]<model_suffix>.bin, use 'x' if no model_suffix present"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[vph,vpv,vsv,vsh,eta].bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. only run in serial mode"
  print '(a)', "  2. vsv = vsv_ref * exp(dlnvs)"
  print '(a)', "  3. elliptic condition: delta = eps, and eta = ((1+eps)*kappa2 - 2)/((1+2*eps)*kappa2 - 2)"
end subroutine


program xsem_model_dlnvs_kappa_thomsen_to_gll

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

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
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vsv_ref
  real(dp), dimension(:,:,:,:), allocatable :: vph,vpv,vsh,vsv,eta
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs,kappa,eps,gamma

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    stop("[ERROR] check your inputs.")
  endif

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

  !====== get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

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

  !====== calculate thomsen parameters
  do iproc = 0, (nproc-1)

    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv_ref', vsv_ref)

    ! read model perturbations
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'dlnvs'//trim(model_suffix), dlnvs)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'kappa'//trim(model_suffix), kappa)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eps'//trim(model_suffix), eps)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'gamma'//trim(model_suffix), gamma)

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        eps(:,:,:,ispec) = 0.0
        gamma(:,:,:,ispec) = 0.0
      endif
    enddo

    vsv = vsv_ref * exp(dlnvs)
    vpv = kappa*vsv
    vph = sqrt((1.0 + 2.0*eps)*vpv**2)
    vsh = sqrt((1.0 + 2.0*gamma)*vsv**2)
    ! elliptic condition: delta = eps
    eta = ((1.0+eps)*kappa**2 - 2.0) / ((1.0+2.0*eps)*kappa**2 - 2.0)

    ! write models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', vph)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', eta)

  enddo

end program
