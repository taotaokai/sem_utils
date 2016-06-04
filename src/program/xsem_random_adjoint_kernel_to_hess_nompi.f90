subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_random_adjoint_kernel_to_hess"
  print '(a)', "    - use random adjoint cijkl/rho kernels to approximate hessian diagonals"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_random_adjoint_kernel_to_hess \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <out_dir>"
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
  print '(a)', "  (string) kernel_dir:  directory holds the kernels files (cijkl_/rho_kernel)"
  print '(a)', "  (string) out_dir:  output directory for hess_kernel"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', "  1. hess_diag ~ sum(K_cijkl,rho**2, for all cijkl,rho)"
  print '(a)', "  2. only run in sequential "

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_random_adjoint_kernel_to_hess

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables
  ! command line args
  integer, parameter :: nargs = 4
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec

  ! kernel gll 
  real(dp), allocatable :: cijkl_kernel(:,:,:,:,:)
  real(dp), allocatable :: rho_kernel(:,:,:,:)
  real(dp), allocatable :: hess_kernel(:,:,:,:)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_reduce_kernel_cijkl_to_lamda_mu: check your inputs."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') kernel_dir
  read(args(4), '(a)') out_dir 

  !====== loop model slices 
  ! get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! initialize gll arrays 
  allocate(cijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(hess_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce cijkl kernels
  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! read kernel gll file
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'cijkl_kernel', cijkl_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'rho_kernel', rho_kernel)

    ! approximate hessian diagonal by summing up squares of all random adjoint kernels
    hess_kernel = sum(cijkl_kernel**2, dim=1) + rho_kernel**2

    print *, "cijkl: min/max=", minval(cijkl_kernel), maxval(cijkl_kernel)
    print *, "rho: min/max=", minval(rho_kernel), maxval(rho_kernel)
    print *, "hess: min/max=", minval(hess_kernel), maxval(hess_kernel)

    ! write out lamda,mu kernel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'hess_kernel', hess_kernel)

  enddo ! iproc

end program
