subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_reduce_kernel_cijkl_to_lamda_mu "
  print '(a)', "    - reduce cijkl kernel to (lamda,mu) kernel"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_cijkl_kernel_to_lamda_mu \"
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
  print '(a)', "  (string) kernel_dir:  directory holds the kernels files"
  print '(a)', "                        proc******_reg1_cijkl_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for lamda,mu_kerenl"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_reduce_kernel_cijkl_to_lamda_mu

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
  real(dp), allocatable :: lamda_kernel(:,:,:,:), mu_kernel(:,:,:,:)

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
  read(args(3), '(a)') kernel_dir
  read(args(4), '(a)') out_dir 

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
  allocate(lamda_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mu_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce cijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read cijkl kernel gll file
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'cijkl_kernel', cijkl_kernel)

    ! reduce cijkl kernel to lamda,mu kernel 
    lamda_kernel = cijkl_kernel(1,:,:,:,:) + &
                   cijkl_kernel(2,:,:,:,:) + &
                   cijkl_kernel(3,:,:,:,:) + &
                   cijkl_kernel(7,:,:,:,:) + &
                   cijkl_kernel(8,:,:,:,:) + &
                   cijkl_kernel(12,:,:,:,:)

    mu_kernel = 2.0*( cijkl_kernel(1,:,:,:,:)  + &
                      cijkl_kernel(7,:,:,:,:)  + &
                      cijkl_kernel(12,:,:,:,:) ) + &
                cijkl_kernel(16,:,:,:,:) + &
                cijkl_kernel(19,:,:,:,:) + &
                cijkl_kernel(21,:,:,:,:)

    print *, "cijkl: min/max=", minval(cijkl_kernel), maxval(cijkl_kernel)
    print *, "lamda: min/max=", minval(lamda_kernel), maxval(lamda_kernel)
    print *, "mu: min/max=", minval(mu_kernel), maxval(mu_kernel)

    ! write out lamda,mu kernel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        'lamda_kernel', lamda_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        'mu_kernel', mu_kernel)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
