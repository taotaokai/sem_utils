subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_cijkl_over_rho_kernel "
  print '(a)', "    - (mijkl, rho_prime) kernel after reparameterization "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_cijkl_kernel_to_lamda_mu \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  mijkl = cijkl / rho"
  print '(a)', "  rho_prime = rho"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir: directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir: directory holds proc*_reg1_vpv,vph,..,rho.bin"
  print '(a)', "  (string) kernel_dir: directory holds proc*_reg1_cijkl,rho_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_cijkl_over_rho_kernel

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
  real(dp), allocatable :: rho_kernel(:,:,:,:)
  ! model gll
  real(dp), dimension(:,:,:,:), allocatable :: vpv, vph, vsv, vsh, eta, rho
  real(dp), dimension(:,:,:,:), allocatable :: A, C, F, L, N

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
  allocate(vpv(NGLLX,NGLLY,NGLLZ,nspec), vph(NGLLX,NGLLY,NGLLZ,nspec), &
           vsv(NGLLX,NGLLY,NGLLZ,nspec), vsh(NGLLX,NGLLY,NGLLZ,nspec), &
           eta(NGLLX,NGLLY,NGLLZ,nspec), rho(NGLLX,NGLLY,NGLLZ,nspec) )
  allocate(A(NGLLX,NGLLY,NGLLZ,nspec), C(NGLLX,NGLLY,NGLLZ,nspec), &
           F(NGLLX,NGLLY,NGLLZ,nspec), L(NGLLX,NGLLY,NGLLZ,nspec), &
           N(NGLLX,NGLLY,NGLLZ,nspec) )
  allocate(cijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce cijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)

    A = vph**2 * rho ! c1, c7
    C = vpv**2 * rho ! c12
    L = vsv**2 * rho ! c16, c19
    N = vsh**2 * rho ! c21
    F = eta * (A - 2.0*L) ! c3, c8

! index of cijkl(1:21,,,,) to Cij(6,6):
!
!   1  2  3| 4  5  6              A A-2N F|
!      7  8| 9 10 11                  A  F|
!        12|13 14 15                     C|
!   -------|--------  => TISO     --------|--------
!          |16 17 18                      |L
!          |   19 20                      |   L
!          |      21                      |      N

    ! read kernels
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'rho_kernel', rho_kernel)
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, cijkl_kernel)

    ! calculate ratio = 1 + (cijkl_kernel*cijkl)/(rho_kernel*rho)
    ! rho_prime_kernel = ratio * rho_kernel
    ! here, re-use eta
    eta = 1.0 + ( cijkl_kernel(1,:,:,:,:)  * A + &
            cijkl_kernel(2,:,:,:,:)  * (A - 2.0*N) + &
            cijkl_kernel(3,:,:,:,:)  * F + &
            cijkl_kernel(7,:,:,:,:)  * A + &
            cijkl_kernel(8,:,:,:,:)  * F + &
            cijkl_kernel(12,:,:,:,:) * C + &
            cijkl_kernel(16,:,:,:,:) * L + &
            cijkl_kernel(19,:,:,:,:) * L + &
            cijkl_kernel(21,:,:,:,:) * N ) / ( rho_kernel * rho )

    ! write out ratio
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho_ratio', eta)
    !call sem_io_write_cijkl_kernel(out_dir, iproc, iregion, 'mijkl' , rho*cijkl_kernel)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
