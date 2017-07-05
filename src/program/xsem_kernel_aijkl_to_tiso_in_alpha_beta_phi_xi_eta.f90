subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_aijkl_to_tiso_in_alpha_beta_phi_xi_eta"
  print '(a)', "    - reduce aijkl kernel to TISO kernel (alpha,beta,phi,xi,eta)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_aijkl_to_tiso_in_alpha_beta_phi_xi_eta \"
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
  print '(a)', "  (string) model_dir:  directory holds the model files, proc*_reg1_[vp0,vs0,alpha,beta,phi,xi,eta].bin"
  print '(a)', "  (string) kernel_dir:  directory holds the aijkl_kernel files"
  print '(a)', "                        proc******_reg1_aijkl_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for [alpha,beta,phi,xi,eta]_kernel"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. aijkl = cijkl/rho"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_aijkl_to_tiso

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
  real(dp), dimension(:,:,:,:), allocatable :: vp0, vs0, alpha, beta
  real(dp), dimension(:,:,:,:), allocatable :: vp, vs, xi, phi, eta

  ! kernel gll 
  real(dp), allocatable :: aijkl_kernel(:,:,:,:,:)
  real(dp), dimension(:,:,:,:), allocatable :: alpha_kernel, beta_kernel
  real(dp), dimension(:,:,:,:), allocatable :: xi_kernel, phi_kernel, eta_kernel

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
  allocate(vp0(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs0(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vp(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(alpha(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(beta(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xi(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(phi(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(aijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))

  allocate(alpha_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(beta_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xi_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(phi_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce aijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read aijkl_kernel
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'aijkl_kernel', aijkl_kernel)

    ! read reference model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vp0', vp0)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vs0', vs0)

    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'alpha', alpha)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'beta', beta)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'phi', phi)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'xi', xi)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)

    ! Voigt average isotropic Vp and Vs
    vp = (1 + alpha)*vp0
    vs = (1 + beta)*vs0

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        xi(:,:,:,ispec) = 0.0
        phi(:,:,:,ispec) = 0.0
        eta(:,:,:,ispec) = 1.0
      endif
    enddo

    ! reduce kernels
    alpha_kernel = ( aijkl_kernel(1,:,:,:,:)*(1.0+1.0/5.0*phi)      &
                   + aijkl_kernel(7,:,:,:,:)*(1.0+1.0/5.0*phi)      &
                   + aijkl_kernel(12,:,:,:,:)*(1.0-4.0/5.0*phi)     &
                   + aijkl_kernel(2,:,:,:,:)*(1.0+1.0/5.0*phi)      &
                   + aijkl_kernel(3,:,:,:,:)*(1.0+1.0/5.0*phi)*eta  &
                   + aijkl_kernel(8,:,:,:,:)*(1.0+1.0/5.0*phi)*eta  &
                   ) * 2.0 * vp0**2 * (1+alpha)

    beta_kernel = ( aijkl_kernel(2,:,:,:,:)*(-2.0*(1.0+2.0/3.0*xi))      &
                  + aijkl_kernel(3,:,:,:,:)*(-2.0*eta*(1.0-1.0/3.0*xi))  &
                  + aijkl_kernel(8,:,:,:,:)*(-2.0*eta*(1.0-1.0/3.0*xi))  &
                  + aijkl_kernel(16,:,:,:,:)*(1.0-1.0/3.0*xi)            &
                  + aijkl_kernel(19,:,:,:,:)*(1.0-1.0/3.0*xi)            &
                  + aijkl_kernel(21,:,:,:,:)*(1.0+2.0/3.0*xi)            &
                  ) * 2.0 * vs0**2 * (1+beta)

    phi_kernel = ( aijkl_kernel(1,:,:,:,:)*1.0/5.0     &
                 + aijkl_kernel(7,:,:,:,:)*1.0/5.0     &
                 + aijkl_kernel(12,:,:,:,:)*(-4.0/5.0) &
                 + aijkl_kernel(2,:,:,:,:)*1.0/5.0     &
                 + aijkl_kernel(3,:,:,:,:)*1.0/5.0*eta &
                 + aijkl_kernel(8,:,:,:,:)*1.0/5.0*eta &
                 ) * vp**2

    xi_kernel = ( aijkl_kernel(16,:,:,:,:)*(-1.0/3.0)        &
                + aijkl_kernel(19,:,:,:,:)*(-1.0/3.0)  &
                + aijkl_kernel(21,:,:,:,:)*(2.0/3.0)  &
                + aijkl_kernel(2,:,:,:,:)*(-4.0/3.0)  &
                + aijkl_kernel(3,:,:,:,:)*(2.0/3.0*eta)  &
                + aijkl_kernel(8,:,:,:,:)*(2.0/3.0*eta)  &
                ) * vs**2

    eta_kernel = ( aijkl_kernel(3,:,:,:,:)*((1.0+1.0/5.0*phi)*vp**2 - 2.0*(1-1.0/3.0*xi)*vs**2) &
                 + aijkl_kernel(8,:,:,:,:)*((1.0+1.0/5.0*phi)*vp**2 - 2.0*(1-1.0/3.0*xi)*vs**2) &
                 )

    ! enforce isotropy for element with ispec_is_tiso = .false.
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        xi_kernel(:,:,:,ispec) = 0.0
        phi_kernel(:,:,:,ispec) = 0.0
        eta_kernel(:,:,:,ispec) = 0.0
      endif
    enddo

    print *, "aijkl_kernel: min/max=", minval(aijkl_kernel), maxval(aijkl_kernel)
    print *, "alpha_kernel: min/max=", minval(alpha_kernel), maxval(alpha_kernel)
    print *, "beta_kernel: min/max=", minval(beta_kernel), maxval(beta_kernel)
    print *, "xi_kernel: min/max=", minval(xi_kernel), maxval(xi_kernel)
    print *, "phi_kernel: min/max=", minval(phi_kernel), maxval(phi_kernel)
    print *, "eta_kernel: min/max=", minval(eta_kernel), maxval(eta_kernel)

    ! write out kernel files
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'alpha_kernel', alpha_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'beta_kernel', beta_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'phi_kernel', phi_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'xi_kernel', xi_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta_kernel', eta_kernel)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
