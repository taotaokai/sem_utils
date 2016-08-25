subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_aijkl_to_vsv2_kappa2_thomsen_elliptic"
  print '(a)', "    - reduce aijkl kernel to VTI kernel (vsv2, kappa2, eps, gamma)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_aijkl_to_vsv2_kappa2_thomsen_elliptic \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  aijkl kernels are indexed from 1 to 21 as defined in "
  print '(a)', "  specfem3d_globe/src/specfem3D/compute_kernels.f90:compute_strain_product()"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vsv2,kappa2,eps,gamma].bin"
  print '(a)', "  (string) kernel_dir:  directory holds the aijkl_kernel files"
  print '(a)', "                        proc******_reg1_aijkl_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for [vsv2,kappa2,eps,gamma]_kernel"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. aijkl = cijkl/rho"
  print '(a)', "  3. thomsen's anisotropic parameters: eps = (C66-C44)/C44/2, gamma = (C11-C33)/C33/2"
  print '(a)', "  4. elliptical condition: delta = eps "

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_aijkl_to_vti_3pars

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
  integer :: nspec, ispec

  ! model gll
  real(dp), dimension(:,:,:,:), allocatable :: vsv2, kappa2, eps, gamma

  ! kernel gll 
  real(dp), allocatable :: aijkl_kernel(:,:,:,:,:)
  real(dp), dimension(:,:,:,:), allocatable :: vsv2_kernel, kappa2_kernel
  real(dp), dimension(:,:,:,:), allocatable :: eps_kernel, gamma_kernel

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_kernel_aijkl_to_vti: check your inputs."
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
  allocate(vsv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappa2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(aijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappa2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eps_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gamma_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce aijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)

    ! read aijkl_kernel
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'aijkl_kernel', aijkl_kernel)

    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv2', vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'kappa2', kappa2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eps', eps)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'gamma', gamma)

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        eps(:,:,:,ispec) = 0.0
        gamma(:,:,:,ispec) = 0.0
      endif
    enddo

    ! reduce kernels
    vsv2_kernel = aijkl_kernel(1,:,:,:,:)*(1.0+2.0*eps)*kappa2                       &
                + aijkl_kernel(7,:,:,:,:)*(1.0+2.0*eps)*kappa2                       &
                + aijkl_kernel(12,:,:,:,:)*kappa2                                    &
                + aijkl_kernel(2,:,:,:,:)*((1.0+2.0*eps)*kappa2-2.0*(1.0+2.0*gamma)) &
                + aijkl_kernel(3,:,:,:,:)*((1.0+eps)*kappa2-2.0)                     &
                + aijkl_kernel(8,:,:,:,:)*((1.0+eps)*kappa2-2.0)                     &
                + aijkl_kernel(16,:,:,:,:)                                           &
                + aijkl_kernel(19,:,:,:,:)                                           &
                + aijkl_kernel(21,:,:,:,:)*(1.0+2.0*gamma)

    kappa2_kernel = aijkl_kernel(1,:,:,:,:)*(1.0+2.0*eps)*vsv2 &
                  + aijkl_kernel(7,:,:,:,:)*(1.0+2.0*eps)*vsv2 &
                  + aijkl_kernel(12,:,:,:,:)*vsv2              &
                  + aijkl_kernel(2,:,:,:,:)*(1.0+2.0*eps)*vsv2 &
                  + aijkl_kernel(3,:,:,:,:)*(1.0+eps)*vsv2     &
                  + aijkl_kernel(8,:,:,:,:)*(1.0+eps)*vsv2                  

    eps_kernel = aijkl_kernel(1,:,:,:,:)*2.0*kappa2*vsv2 &
               + aijkl_kernel(7,:,:,:,:)*2.0*kappa2*vsv2 &
               + aijkl_kernel(2,:,:,:,:)*2.0*kappa2*vsv2 &
               + aijkl_kernel(3,:,:,:,:)*kappa2*vsv2     &
               + aijkl_kernel(8,:,:,:,:)*kappa2*vsv2

    gamma_kernel = aijkl_kernel(21,:,:,:,:)*2.0*vsv2 &
                 - aijkl_kernel(2,:,:,:,:)*4.0*vsv2

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        eps_kernel(:,:,:,ispec) = 0.0
        gamma_kernel(:,:,:,ispec) = 0.0
      endif
    enddo

    print *, "aijkl_kernel: min/max=", minval(aijkl_kernel), maxval(aijkl_kernel)
    print *, "vsv2_kernel: min/max=", minval(vsv2_kernel), maxval(vsv2_kernel)
    print *, "kappa2_kernel: min/max=", minval(kappa2_kernel), maxval(kappa2_kernel)
    print *, "eps_kernel: min/max=", minval(eps_kernel), maxval(eps_kernel)
    print *, "gamma_kernel: min/max=", minval(gamma_kernel), maxval(gamma_kernel)

    ! write out kernel files
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv2_kernel', vsv2_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'kappa2_kernel', kappa2_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eps_kernel', eps_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'gamma_kernel', gamma_kernel)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
