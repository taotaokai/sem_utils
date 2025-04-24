subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_cov_event_kernels_1"
  print '(a)', "    - calculate covariance matrix of event kernels"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_cov_event_kernels_1 \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir_list> <kernel_name> "
  print '(a)', "    <out_file>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir_list:  list of directories containing kernel files"
  print '(a)', "  (string) kernel_tags:  e.g. vp2_kernel,vsv2_kernel,vsh2_kernel"
  print '(a)', "  (string) out_file:  output file for covariance matrix"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', ""

end subroutine


!//////////////////////////////////////////////////////////////////////////////
program xsem_cov_event_kernels_1

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
  character(len=MAX_STRING_LEN) :: kernel_dir_list
  character(len=MAX_STRING_LEN) :: kernel_tags
  character(len=MAX_STRING_LEN) :: out_file

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank
  ! list of kernel directories
  character(len=MAX_STRING_LEN), allocatable :: kernel_dirs(:)
  integer :: iker, nkernel
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! kernel gll 
  real(dp), allocatable :: kernel_gll(:,:,:,:,:)
  ! covariance matrix
  real(dp), allocatable :: cov(:,:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_cov_event_kernels_1: check your inputs."
      call abort_mpi()
    endif 
  endif
  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') kernel_dir_list
  read(args(4), '(a)') kernel_name
  read(args(5), '(a)') out_file

  call synchronize_all()

  !====== read kernel_dir_list
  call sem_utils_read_line(kernel_dir_list, kernel_dirs, nkernel)

  !====== loop model slices 

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize gll arrays
  allocate(kernel(nkernel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(cov(nkernel,nkernel))

  ! combine event kernels
  do iproc = myrank, (nproc-1), nrank

    print *, '#-- iproc=', iproc

    kernel_sum = 0.0_dp
    do iker = 1, nkernel

      ! read kernel gll
      call sem_io_read_gll_file_n(kernel_dirs(iker), iproc, iregion, &
        kernel_name, kernel)

      ! read mask
      if (use_mask) then
        call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, &
          mask_tag, mask)
      endif

      ! sum up kernels
      kernel_sum = kernel_sum + kernel * mask

    enddo ! iker

    ! write out gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
      kernel_name, kernel_sum)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
