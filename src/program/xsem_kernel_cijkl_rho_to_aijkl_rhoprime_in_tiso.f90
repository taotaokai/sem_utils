subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso "
  print '(a)', "    - transform from cijkl,rho to aijkl,rhoprime kernel in TISO model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  cijkl kernels are indexed from 1 to 21 as defined in "
  print '(a)', "  specfem3d_globe/src/specfem3D/compute_kernels.f90:compute_strain_product()"
  print '(a)', ""
  print '(a)', "  aijkl = cijkl/rho; rhoprime = rho"
  print '(a)', ""
  print '(a)', "  aijkl_kernel = rho * cijkl_kernel"
  print '(a)', "  rhoprime_kernel = rho_kernel + sum(cijkl * cijkl_kernel, over ijkl)/rho"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir: proc*_reg1_[vpv|vph|vsv|vsh|eta|rho].bin"
  print '(a)', "  (string) kernel_dir: proc*_reg1_[cijkl|rho]_kernel.bin"
  print '(a)', "  (string) out_dir: proc*_reg1_[aijkl|rhoprime]_kernel.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_cijkl_rho_to_aijkl_rhoprime_in_tiso

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

  ! model gll
  real(dp), allocatable :: cijkl(:,:,:,:,:)
  real(dp), dimension(:,:,:,:), allocatable :: rho, gll1

  ! kernel gll 
  real(dp), allocatable :: cijkl_kernel(:,:,:,:,:)
  real(dp), allocatable :: rho_kernel(:,:,:,:)

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
  allocate(cijkl(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll1(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(cijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce cijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    !--- get density model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)

    !--- get cijkl model, converted from tiso model files 
    cijkl = 0.0
    ! A = vph**2 * rho; 
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', gll1)
    cijkl(1,:,:,:,:) = gll1**2 * rho
    cijkl(7,:,:,:,:) = cijkl(1,:,:,:,:)
    ! C = vpv**2 * rho;
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', gll1)
    cijkl(12,:,:,:,:) = gll1**2 * rho
    ! L = vsv**2 * rho; 
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', gll1)
    cijkl(16,:,:,:,:) = gll1**2 * rho
    cijkl(19,:,:,:,:) = cijkl(16,:,:,:,:)
    ! N = vsh**2 * rho;
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', gll1)
    cijkl(21,:,:,:,:) = gll1**2 * rho
    ! F = eta * (A - 2.0*L)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', gll1)
    cijkl(3,:,:,:,:) = gll1 * (cijkl(1,:,:,:,:) - 2.0*cijkl(16,:,:,:,:))
    cijkl(8,:,:,:,:) = cijkl(3,:,:,:,:)
    ! cijkl(2,::::) = A - 2.0*N
    cijkl(2,:,:,:,:) = cijkl(1,:,:,:,:) - 2.0*cijkl(21,:,:,:,:)

    !--- get cijkl,rho kernel
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'rho_kernel', rho_kernel)
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'cijkl_kernel', cijkl_kernel)

    !--- transform to aijkl,rhoprime kernel
    ! aijkl_kernel = cijkl_kernel * rho
    cijkl_kernel = cijkl_kernel * spread(rho, 1, 21)
    ! rhoprime_kernel = rho_kernel + sum(aijkl_kernel*cijkl, over ijkl)/rho**2
    ! note: in SEM the kernel includes the factor in reduction from ijkl to 21 indices  
    rho_kernel = rho_kernel + sum(cijkl*cijkl_kernel, dim=1)/rho**2

    print *, "aijkl: min/max=", minval(cijkl_kernel), maxval(cijkl_kernel)
    print *, "rhoprime: min/max=", minval(rho_kernel), maxval(rho_kernel)

    !--- write out aijkl,rhoprime kernel
    call sem_io_write_cijkl_kernel(out_dir, iproc, iregion, &
        'aijkl_kernel', cijkl_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        'rhoprime_kernel', rho_kernel)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
