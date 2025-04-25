subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_alpha_beta_phi_xi_eta"
  print '(a)', "    - add model perturbation "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_alpha_beta_phi_xi_eta"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <dmodel_dir> <dmodel_suffix>"
  print '(a)', "    <step_length> "
  print '(a)', "    <min_alpha> <max_alpha>"
  print '(a)', "    <min_beta> <max_beta>"
  print '(a)', "    <min_phi> <max_phi>"
  print '(a)', "    <min_xi> <max_xi>"
  print '(a)', "    <min_eta> <max_eta>"
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
  print '(a)', "  (string) model_dir:  directory holds old model files proc000***_reg1_[alpha,beta,phi,xi,eta].bin"
  print '(a)', "  (string) dmodel_dir:  directory holds dmodel files"
  print '(a)', "  (string) dmodel_suffix:  proc000***_reg1_[alpha,beta,phi,xi,eta]<dmodel_suffix>.bin"
  print '(a)', "  (float) step_length:  step size multiplied on dmodel"
  print '(a)', "  (float) min_alpha:  minimum value"
  print '(a)', "  (float) max_alpha:  maximum value"
  print '(a)', "  (float) min_beta:  minimum value"
  print '(a)', "  (float) max_beta:  maximum value"
  print '(a)', "  (float) min_phi:    minimum value"
  print '(a)', "  (float) max_phi:    maximum value"
  print '(a)', "  (float) min_xi:  minimum value"
  print '(a)', "  (float) max_xi:  maximum value"
  print '(a)', "  (float) min_eta:  minimum value"
  print '(a)', "  (float) max_eta:  maximum value"
  print '(a)', "  (int) output_dmodel:  flag out dmodel (0:no, 1:yes) proc000***_reg1_[alpha,beta,phi,xi,eta]_dmodel.bin"
  print '(a)', "  (string) out_dir:  out directory for new model proc000***_reg1_[alpha,beta,phi,xi,eta].bin"
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
  integer, parameter :: nargs = 18
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: dmodel_dir
  character(len=MAX_STRING_LEN) :: dmodel_suffix
  real(dp) :: step_length
  real(dp) :: min_alpha, min_beta, min_phi, min_xi, min_eta
  real(dp) :: max_alpha, max_beta, max_phi, max_xi, max_eta
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
  real(dp), dimension(:,:,:,:), allocatable :: alpha, beta, phi, xi, eta
  real(dp), dimension(:,:,:,:), allocatable :: alpha_new, beta_new, phi_new, xi_new, eta_new
  ! dmodel 
  real(dp), dimension(:,:,:,:), allocatable :: alpha_dmodel, beta_dmodel, phi_dmodel, xi_dmodel, eta_dmodel 

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
  read(args(7),*) min_alpha
  read(args(8),*) max_alpha
  read(args(9),*) min_beta
  read(args(10),*) max_beta
  read(args(11),*) min_phi
  read(args(12),*) max_phi
  read(args(13),*) min_xi
  read(args(14),*) max_xi 
  read(args(15),*) min_eta
  read(args(16),*) max_eta 
  read(args(17),*) output_dmodel
  read(args(18),'(a)') out_dir

  ! process input
  if (dmodel_suffix == 'x') then
    dmodel_suffix = ''
  endif
  if ( min_beta>max_beta .or. min_alpha>max_alpha &
      .or. min_phi>max_phi .or. min_xi>max_xi .or. min_eta>max_eta) then
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
  allocate(alpha(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(beta(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(phi(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xi(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(alpha_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(beta_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(phi_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xi_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta_new(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(alpha_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(beta_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(phi_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xi_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== create new model
  do iproc = myrank, (nproc-1), nrank

    print *, "====== proc ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read old models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'alpha', alpha)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'beta', beta)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'phi', phi)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'xi', xi)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)

    ! read dmodel
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'alpha'//trim(dmodel_suffix), alpha_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'beta'//trim(dmodel_suffix), beta_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'phi'//trim(dmodel_suffix), phi_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'xi'//trim(dmodel_suffix), xi_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'eta'//trim(dmodel_suffix), eta_dmodel)

    ! add dmodel with step_length
    alpha_new = alpha + step_length*alpha_dmodel
    beta_new = beta + step_length*beta_dmodel
    phi_new = phi   + step_length*phi_dmodel
    xi_new = xi + step_length*xi_dmodel
    eta_new = eta + step_length*eta_dmodel

    ! limit model value range
    where(alpha_new < min_alpha) alpha_new = min_alpha
    where(alpha_new > max_alpha) alpha_new = max_alpha
    where(beta_new < min_beta)   beta_new = min_beta
    where(beta_new > max_beta)   beta_new = max_beta
    where(phi_new < min_phi)     phi_new = min_phi
    where(phi_new > max_phi)     phi_new = max_phi
    where(xi_new < min_xi)       xi_new = min_xi
    where(xi_new > max_xi)       xi_new = max_xi
    where(eta_new < min_eta)     eta_new = min_eta
    where(eta_new > max_eta)     eta_new = max_eta

    ! enforce isotropic element
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        phi(:,:,:,ispec) = 0.0
        xi(:,:,:,ispec) = 0.0
        eta(:,:,:,ispec) = 1.0

        phi_new(:,:,:,ispec) = 0.0
        xi_new(:,:,:,ispec) = 0.0
        eta_new(:,:,:,ispec) = 1.0
      endif
    enddo

    print *, "alpha min/max = ", minval(alpha_new), maxval(alpha_new)
    print *, "beta min/max = ", minval(beta_new), maxval(beta_new)
    print *, "phi min/max = ",   minval(phi_new),   maxval(phi_new)
    print *, "xi min/max = ", minval(xi_new), maxval(xi_new)
    print *, "eta min/max = ", minval(eta_new), maxval(eta_new)

    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'alpha', alpha_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'beta', beta_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'phi', phi_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'xi', xi_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', eta_new)

    ! write out dmodel
    if (output_dmodel == 1) then
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'alpha_dmodel', alpha_new - alpha)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'beta_dmodel', beta_new - beta)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'phi_dmodel',   phi_new -   phi)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'xi_dmodel', xi_new - xi)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta_dmodel', eta_new - eta)
    endif

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
