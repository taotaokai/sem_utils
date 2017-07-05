subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_scale_dmodel_alpha_beta_phi_xi_eta"
  print '(a)', "    - scale model perturbation "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_scale_dmodel "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <model_dir> <model_suffix> "
  print '(a)', "    <max_beta> "
  print '(a)', "    <out_dir> <out_suffix>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  max(abs(dmodel_scaled/model)) = max_beta"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds model files"
  print '(a)', "  (string) model_suffix:  proc000***_reg1_[alpha,beta,phi,xi,eta]<model_suffix>.bin"
  print '(a)', "  (float) max_beta:  maximum abs(beta) "
  print '(a)', "  (string) out_dir:  out directory "
  print '(a)', "  (string) out_suffix:  proc000***_reg1_[alpha,beta,phi,xi,eta]<out_suffix>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. scaling factor is determined to have the maximum amplitude of scaled beta dmodel <max_beta>"
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
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_suffix
  real(dp) :: max_beta
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_suffix

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! dmodel
  real(dp), dimension(:,:,:,:), allocatable :: dmodel
  ! scale factor
  real(dp) :: max_beta_local, max_beta_all
  real(dp) :: scale_factor

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
  read(args(5),*) max_beta
  read(args(6),'(a)') out_dir
  read(args(7),'(a)') out_suffix

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(dmodel(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== get scale factor
  max_beta_local = 0.0
  do iproc = myrank, (nproc-1), nrank
    ! read dmodel
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, "beta"//trim(model_suffix), dmodel)
    ! get max amplitude 
    max_beta_local = max(max_beta_local, maxval(abs(dmodel)))
  enddo
  call max_all_dp(max_beta_local, max_beta_all)

  if (myrank == 0) then
    scale_factor = max_beta/max_beta_all
  endif
  call bcast_all_singledp(scale_factor)

  print *, "scale_factor = ", scale_factor

  !====== scale dmodel
  do iproc = myrank, (nproc-1), nrank

    print *, "======", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! beta
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, "beta"//trim(model_suffix), dmodel)
    dmodel = scale_factor*dmodel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, "beta"//trim(out_suffix), dmodel)
    print *, "beta min/max = ", minval(dmodel), maxval(dmodel)

    ! alpha 
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, "alpha"//trim(model_suffix), dmodel)
    dmodel = scale_factor*dmodel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, "alpha"//trim(out_suffix), dmodel)
    print *, "alpha min/max = ", minval(dmodel), maxval(dmodel)

    ! phi
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, "phi"//trim(model_suffix), dmodel)
    dmodel = scale_factor*dmodel
    ! enforce isotropic elements
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        dmodel(:,:,:,ispec) = 0.0
      endif
    enddo
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, "phi"//trim(out_suffix), dmodel)
    print *, "phi min/max = ", minval(dmodel), maxval(dmodel)

    ! xi
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, "xi"//trim(model_suffix), dmodel)
    dmodel = scale_factor*dmodel
    ! enforce isotropic elements
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        dmodel(:,:,:,ispec) = 0.0
      endif
    enddo
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, "xi"//trim(out_suffix), dmodel)
    print *, "xi min/max = ", minval(dmodel), maxval(dmodel)

    ! eta
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, "eta"//trim(model_suffix), dmodel)
    dmodel = scale_factor*dmodel
    ! enforce isotropic elements
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        dmodel(:,:,:,ispec) = 0.0
      endif
    enddo
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, "eta"//trim(out_suffix), dmodel)
    print *, "eta min/max = ", minval(dmodel), maxval(dmodel)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
