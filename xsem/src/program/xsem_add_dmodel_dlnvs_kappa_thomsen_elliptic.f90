subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_dlnvs_kappa_thomsen_elliptic"
  print '(a)', "    - add model perturbation "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_dlnvs_kappa_thomsen_elliptic"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <dmodel_dir> <dmodel_suffix>"
  print '(a)', "    <step_length> "
  print '(a)', "    <min_dlnvs> <max_dlnvs>"
  print '(a)', "    <min_kappa> <max_kappa>"
  print '(a)', "    <min_eps> <max_eps>"
  print '(a)', "    <min_gamma> <max_gamma>"
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
  print '(a)', "  (string) model_dir:  directory holds old model files proc000***_reg1_[dlnvs,kappa,eps,gamma].bin"
  print '(a)', "  (string) dmodel_dir:  directory holds dmodel files"
  print '(a)', "  (string) dmodel_suffix:  proc000***_reg1_[dlnvs,kappa,eps,gamma]<dmodel_suffix>.bin"
  print '(a)', "  (float) step_length:  step size multiplied on dmodel"
  print '(a)', "  (float) min_dlnvs:  minimum value"
  print '(a)', "  (float) max_dlnvs:  maximum value"
  print '(a)', "  (float) min_kappa:  minimum value"
  print '(a)', "  (float) max_kappa:  maximum value"
  print '(a)', "  (float) min_eps:    minimum value"
  print '(a)', "  (float) max_eps:    maximum value"
  print '(a)', "  (float) min_gamma:  minimum value"
  print '(a)', "  (float) max_gamma:  maximum value"
  print '(a)', "  (int) output_dmodel:  flag out dmodel (0:no, 1:yes) proc000***_reg1_[dlnvs,kappa,eps,gamma]_dmodel.bin"
  print '(a)', "  (string) out_dir:  out directory for new model proc000***_reg1_[dlnvs,kappa,eps,gamma].bin"
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
  integer, parameter :: nargs = 16
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: dmodel_dir
  character(len=MAX_STRING_LEN) :: dmodel_suffix
  real(dp) :: step_length
  real(dp) :: min_dlnvs, min_kappa, min_eps, min_gamma
  real(dp) :: max_dlnvs, max_kappa, max_eps, max_gamma
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
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs, kappa, eps, gamma
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs_new, kappa_new, eps_new, gamma_new
  ! dmodel 
  real(dp), dimension(:,:,:,:), allocatable :: dlnvs_dmodel, kappa_dmodel, eps_dmodel, gamma_dmodel 

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
  read(args(7),*) min_dlnvs
  read(args(8),*) max_dlnvs
  read(args(9),*) min_kappa
  read(args(10),*) max_kappa
  read(args(11),*) min_eps
  read(args(12),*) max_eps
  read(args(13),*) min_gamma
  read(args(14),*) max_gamma 
  read(args(15),*) output_dmodel
  read(args(16),'(a)') out_dir

  ! process input
  if (dmodel_suffix == 'x') then
    dmodel_suffix = ''
  endif
  if ( min_dlnvs>max_dlnvs .or. min_kappa>max_kappa &
      .or. min_eps>max_eps .or. min_gamma>max_gamma) then
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
  allocate(dlnvs(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappa(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dlnvs_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappa_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(  eps_new(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma_new(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dlnvs_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappa_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma_dmodel(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== create new model
  do iproc = myrank, (nproc-1), nrank

    print *, "====== proc ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read old models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'dlnvs', dlnvs)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'kappa', kappa)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eps', eps)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'gamma', gamma)

    ! read dmodel
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'dlnvs'//trim(dmodel_suffix), dlnvs_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'kappa'//trim(dmodel_suffix), kappa_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion,   'eps'//trim(dmodel_suffix), eps_dmodel)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'gamma'//trim(dmodel_suffix), gamma_dmodel)

    ! add dmodel with step_length
    dlnvs_new = dlnvs + step_length*dlnvs_dmodel
    kappa_new = kappa + step_length*kappa_dmodel
      eps_new = eps   + step_length*eps_dmodel
    gamma_new = gamma + step_length*gamma_dmodel

    ! limit model value range
    where(dlnvs_new < min_dlnvs) dlnvs_new = min_dlnvs
    where(dlnvs_new > max_dlnvs) dlnvs_new = max_dlnvs
    where(kappa_new < min_kappa) kappa_new = min_kappa
    where(kappa_new > max_kappa) kappa_new = max_kappa
    where(  eps_new < min_eps)     eps_new = min_eps
    where(  eps_new > max_eps)     eps_new = max_eps
    where(gamma_new < min_gamma) gamma_new = min_gamma
    where(gamma_new > max_gamma) gamma_new = max_gamma

    ! enforce isotropic element
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        eps(:,:,:,ispec) = 0.0
        gamma(:,:,:,ispec) = 0.0
        eps_new(:,:,:,ispec) = 0.0
        gamma_new(:,:,:,ispec) = 0.0
      endif
    enddo

    print *, "dlnvs min/max = ", minval(dlnvs_new), maxval(dlnvs_new)
    print *, "kappa min/max = ", minval(kappa_new), maxval(kappa_new)
    print *, "eps min/max = ",   minval(eps_new),   maxval(eps_new)
    print *, "gamma min/max = ", minval(gamma_new), maxval(gamma_new)

    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dlnvs', dlnvs_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'kappa', kappa_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion,   'eps',   eps_new)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'gamma', gamma_new)

    ! write out dmodel
    if (output_dmodel == 1) then
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'dlnvs_dmodel', dlnvs_new - dlnvs)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'kappa_dmodel', kappa_new - kappa)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion,   'eps_dmodel',   eps_new -   eps)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'gamma_dmodel', gamma_new - gamma)
    endif

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
