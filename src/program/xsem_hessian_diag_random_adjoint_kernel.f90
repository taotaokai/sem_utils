subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_hessian_diag_random_adjoint_kernel"
  print '(a)', "    - use random adjoint cijkl kernels to approximate hessian diagonals"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_random_adjoint_kernel_to_hess \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  cijkl kernels are indexed from 1 to 21 as defined in "
  print '(a)', "  specfem3d_globe/src/specfem3D/compute_kernels.f90:compute_strain_product()"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_rho.bin"
  print '(a)', "  (string) kernel_dir:  directory holds the kernels files (cijkl_kernel)"
  print '(a)', "  (string) out_dir:  output directory for hess_kernel"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. for Aijkl hessian: hess_diag ~ rho^2*sum(K_cijkl**2), ijkl is chosen for P, S and PS, see code"


end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_random_adjoint_kernel_to_hess

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
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! mpi
  integer :: myrank, nrank

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec

  ! kernel gll 
  real(dp), allocatable :: cijkl_kernel(:,:,:,:,:)
  real(dp), allocatable :: hess(:,:,:,:)
  !real(dp), allocatable :: hess_p(:,:,:,:)
  !real(dp), allocatable :: hess_s(:,:,:,:)
  !real(dp), allocatable :: hess_ps(:,:,:,:)

  ! model gll
  real(dp), allocatable :: rho(:,:,:,:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_reduce_kernel_cijkl_to_lamda_mu: check your inputs."
      call abort_mpi()
    endif 
  endif
  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') kernel_dir
  read(args(5), '(a)') out_dir 

  !====== loop model slices 

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize gll arrays 
  allocate(cijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(hess(NGLLX,NGLLY,NGLLZ,nspec))
  !allocate(hess_p(NGLLX,NGLLY,NGLLZ,nspec))
  !allocate(hess_s(NGLLX,NGLLY,NGLLZ,nspec))
  !allocate(hess_ps(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce cijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read kernel gll file
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'cijkl_kernel', cijkl_kernel)
    !call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'rho_kernel', rho_kernel)

    ! read rho gll
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)

    ! approximate hessian diagonal by summing up squares of all random adjoint kernels
    !hess_p = cijkl_kernel(1,:,:,:,:)**2  &
    !       + cijkl_kernel(2,:,:,:,:)**2  &
    !       + cijkl_kernel(3,:,:,:,:)**2  &
    !       + cijkl_kernel(7,:,:,:,:)**2  &
    !       + cijkl_kernel(8,:,:,:,:)**2  &
    !       + cijkl_kernel(12,:,:,:,:)**2

    !hess_s = cijkl_kernel(16,:,:,:,:)**2 &
    !       + cijkl_kernel(17,:,:,:,:)**2 &
    !       + cijkl_kernel(18,:,:,:,:)**2 &
    !       + cijkl_kernel(19,:,:,:,:)**2 &
    !       + cijkl_kernel(20,:,:,:,:)**2 &
    !       + cijkl_kernel(21,:,:,:,:)**2

    !hess_ps = cijkl_kernel(4,:,:,:,:)**2  &
    !        + cijkl_kernel(5,:,:,:,:)**2  &
    !        + cijkl_kernel(6,:,:,:,:)**2  &
    !        + cijkl_kernel(9,:,:,:,:)**2  &
    !        + cijkl_kernel(10,:,:,:,:)**2 &
    !        + cijkl_kernel(11,:,:,:,:)**2 &
    !        + cijkl_kernel(13,:,:,:,:)**2 &
    !        + cijkl_kernel(14,:,:,:,:)**2 &
    !        + cijkl_kernel(15,:,:,:,:)**2

    hess = sum(cijkl_kernel(:,:,:,:,:)**2, dim=1) * rho**2

    print *, "cijkl: min/max=", minval(cijkl_kernel), maxval(cijkl_kernel)
    print *, "hess: min/max=", minval(hess), maxval(hess)
    !print *, "hess_p: min/max=", minval(hess_p), maxval(hess_p)
    !print *, "hess_s: min/max=", minval(hess_s), maxval(hess_s)
    !print *, "hess_ps: min/max=", minval(hess_ps), maxval(hess_ps)

    ! write out diagonal hessian 
    !call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'hess_p', hess_p)
    !call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'hess_s', hess_s)
    !call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'hess_ps', hess_ps)
    ! sum of all hess
    !call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'hess', hess_p+hess_s+hess_ps)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'hess', hess)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
