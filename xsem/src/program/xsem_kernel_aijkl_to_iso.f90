subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_aijkl_to_iso"
  print '(a)', "    - reduce aijkl kernel to isotropic (vp2, vs2) kernel"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_aijkl_to_iso \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <out_dir>"
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
  print '(a)', "  (string) kernel_dir:  directory holds the aijkl_kernel files"
  print '(a)', "                        proc******_reg1_aijkl_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for vp2,vs2_kernel"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. vp2 = vp**2 = (lamda + 2*mu)/rho, vs2 = vs**2 = mu/rho"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_aijkl_to_iso

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
  real(dp), allocatable :: aijkl_kernel(:,:,:,:,:)
  real(dp), allocatable :: vp2_kernel(:,:,:,:), vs2_kernel(:,:,:,:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_kernel_aijkl_to_iso: check your inputs."
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
  allocate(aijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vp2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs2_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! reduce cijkl kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read cijkl kernel gll file
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, 'aijkl_kernel', aijkl_kernel)

    ! reduce cijkl kernel to a,b kernel 
    vp2_kernel = aijkl_kernel(1,:,:,:,:)  &
               + aijkl_kernel(7,:,:,:,:)  &
               + aijkl_kernel(12,:,:,:,:) &
               + aijkl_kernel(2,:,:,:,:)  &
               + aijkl_kernel(3,:,:,:,:)  &
               + aijkl_kernel(8,:,:,:,:)

    vs2_kernel = -2.d0*aijkl_kernel(2,:,:,:,:) &
               - 2.d0*aijkl_kernel(3,:,:,:,:)  &
               - 2.d0*aijkl_kernel(8,:,:,:,:)  &
               + aijkl_kernel(16,:,:,:,:)      &
               + aijkl_kernel(19,:,:,:,:)      &
               + aijkl_kernel(21,:,:,:,:)

    print *, "aijkl_kernel: min/max=", minval(aijkl_kernel), maxval(aijkl_kernel)
    print *, "vp2_kernel: min/max=", minval(vp2_kernel), maxval(vp2_kernel)
    print *, "vs2_kernel: min/max=", minval(vs2_kernel), maxval(vs2_kernel)

    ! write out a,b kernel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vp2_kernel', vp2_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vs2_kernel', vs2_kernel)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
