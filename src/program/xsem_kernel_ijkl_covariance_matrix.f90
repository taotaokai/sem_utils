subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_ijkl_convariance_matrix"
  print '(a)', "    - get covariance matrix of cijkl/aijkl kernel for principal component analysis"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_ijkl_convariance_matrix \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <kernel_tag> <out_file>"
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
  print '(a)', "                        proc******_reg1_<kernel_tag>.bin"
  print '(a)', "  (string) kernel_tag:  kernel file tags"
  print '(a)', "  (string) out_file:  output file for covariance matrix (21, 21)"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_ijkl_convariance_matrix

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
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: kernel_tag
  character(len=MAX_STRING_LEN) :: out_file

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, j, ier, iproc

  ! mpi
  integer :: myrank, nrank

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec

  ! gll 
  real(dp), allocatable :: kernel_ijkl(:,:,:,:,:)
  real(dp), allocatable :: kernel_cov(:,:), kernel_cov_sum(:,:)
  real(dp), allocatable :: kernel_mean(:), kernel_mean_sum(:)
  real(dp), allocatable :: volume_gll(:,:,:,:)
  real(dp) :: volume, volume_sum

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
  read(args(3), '(a)') kernel_dir
  read(args(4), '(a)') kernel_tag
  read(args(5), '(a)') out_file

  !====== loop model slices 

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize gll arrays 
  allocate(kernel_ijkl(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kernel_cov(21,21), kernel_cov_sum(21,21))
  allocate(kernel_mean(21), kernel_mean_sum(21))
  allocate(volume_gll(NGLLX,NGLLY,NGLLZ,nspec))

  ! calculate kernel mean, variance
  volume = 0.0_dp
  kernel_mean = 0.0_dp
  kernel_cov = 0.0_dp

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! get integral volumes for each gll points
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    call sem_mesh_gll_volume(mesh_data, volume_gll)

    ! read kernel file
    call sem_io_read_cijkl_kernel(kernel_dir, iproc, iregion, kernel_tag, kernel_ijkl)

    ! calculate mean (no need to divide the volume)
    do i = 1, 21
      kernel_mean(i) = kernel_mean(i) + sum(kernel_ijkl(i,:,:,:,:) * volume_gll)
    enddo

    ! calculate variance
    do i = 1, 21
      do j = 1, 21
        kernel_cov(i,j) = kernel_cov(i,j) + &
          sum(kernel_ijkl(i,:,:,:,:) * kernel_ijkl(j,:,:,:,:) * volume_gll)
      enddo
    enddo

    ! local volume
    volume = volume + sum(volume_gll)

  enddo ! iproc

  ! sum up over all processes 
  call synchronize_all()
  call sum_all_dp(volume, volume_sum)
  call sum_all_dp_n(kernel_mean, kernel_mean_sum, 21)
  call sum_all_dp_n(kernel_cov, kernel_cov_sum, 21**2)

  ! write out kernel_mean,cov
  if (myrank == 0) then
    open(IOUT, file=trim(out_file), status='unknown', form='formatted', action='write', iostat=ier)

    if (ier /= 0) then
      write(*,*) '[ERROR] failed to open file ', trim(out_file)
      call abort_mpi()
    endif

    kernel_mean_sum = kernel_mean_sum/volume_sum
    write(IOUT, '(a)') "#kernel_mean"
    write(IOUT, '( 21(ES15.7,2X) )') kernel_mean_sum 

    kernel_cov_sum = kernel_cov_sum/volume_sum
    do i = 1, 21
      do j = 1, 21
        kernel_cov_sum(i,j) = kernel_cov_sum(i,j) - kernel_mean_sum(i)*kernel_mean_sum(j)
      enddo
    enddo

    write(IOUT, '(a)') "#kernel_cov"
    do i = 1, 21
      write(IOUT, '( 21(ES15.7,2X) )') kernel_cov_sum(i, :)
    enddo

  endif 

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
