subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_1"
  print '(a)', "    - sum up event kernels"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_1 \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir_list> <kernel_name> <use_mask>"
  print '(a)', "    <nroot> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir_list:  list of event kernel directories (proc***_reg1_***_kernel.bin)"
  print '(a)', "  (string) kernel_name:  kernel name to be summed (e.g. rho_kernel)"
  print '(a)', "  (int) use_mask:  flag whether apply kernel mask (mask file is read from each kernel_dir)"
  print '(a)', "  (int) nroot:  n-th root stacking (must be positive integer)"
  print '(a)', "  (string) out_dir:  directory to write out summed kernel files"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


!//////////////////////////////////////////////////////////////////////////////
program xsem_sum_event_kernels_cijkl

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir_list
  character(len=MAX_STRING_LEN) :: kernel_name
  logical :: use_mask
  integer :: nroot
  character(len=MAX_STRING_LEN) :: out_dir

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
  ! mask
  real(dp), allocatable :: mask(:,:,:,:)
  ! kernel gll 
  real(dp), allocatable :: kernel(:,:,:,:)
  real(dp), allocatable :: kernel_sum(:,:,:,:)
  ! n-th root stacking
  real(dp) :: one_over_nroot
  real(dp), allocatable :: sign_kernel(:,:,:,:)
  real(dp), allocatable :: ones_kernel(:,:,:,:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_sum_event_kernels_1: check your inputs."
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
  select case (args(5))
    case ('0')
      use_mask = .false.
    case ('1')
      use_mask = .true.
    case default
      if (myrank==0) then
        print *, '[ERROR]: use_mask must be 0 or 1'
        call abort_mpi()
      endif
  end select
  read(args(6), *) nroot
  if (nroot < 1) then
    if (myrank==0) then
      print *, '[ERROR]: <nroot> must g.e. 1'
      call abort_mpi()
    endif
  endif
  read(args(7), '(a)') out_dir 

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
  allocate(kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kernel_sum(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mask(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(sign_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ones_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! combine event kernels
  one_over_nroot = 1.0_dp / nroot
  mask = 1.0_dp
  sign_kernel = 1.0_dp
  ones_kernel = 1.0_dp
  do iproc = myrank, (nproc-1), nrank

    print *, '#-- iproc=', iproc

    kernel_sum = 0.0_dp
    do iker = 1, nkernel

      ! read kernel gll
      call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, &
        kernel_name, kernel)

      ! read mask
      if (use_mask) then
        call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, &
          "mask", mask)
      endif

      ! n-th root stacking
      sign_kernel = sign(ones_kernel, kernel)
      kernel_sum = kernel_sum &
        + sign_kernel * abs(kernel)**one_over_nroot * mask

    enddo ! iker

    ! n-th root stacking
    sign_kernel = sign(ones_kernel, kernel_sum)
    kernel_sum = sign_kernel * abs(kernel_sum)**nroot

    ! write out gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
      kernel_name, kernel_sum)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
