subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_cijkl"
  print '(a)', "    - sum up event kernels"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_cijkl \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir_list> <use_mask> "
  print '(a)', "    <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  "
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir_list:  list of event kernel directories (proc*_reg1_cijkl_kernel.bin)"
  print '(a)', "  (int) use_mask:  flag whether apply kernel mask (mask file is read from each kernel_dir)"
  print '(a)', "  (string) mask_tag:  tags of mask files"
  print '(a)', "  (string) out_dir:  output directory of summed cijkl kernel"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. mask files: <kernel_dir>/proc***_reg1_<mask_tag>.bin"
  print '(a)', ""

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
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir_list
  logical :: use_mask
  character(len=MAX_STRING_LEN) :: mask_tag 
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
  real(dp), allocatable :: cijkl_kernel(:,:,:,:,:)
  real(dp), allocatable :: cijkl_kernel_sum(:,:,:,:,:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_sum_event_kernels_cijkl: check your inputs."
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
  select case (args(4))
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
  read(args(5), '(a)') mask_tag 
  read(args(6), '(a)') out_dir

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
  allocate(cijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(cijkl_kernel_sum(21,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mask(NGLLX,NGLLY,NGLLZ,nspec))

  ! combine event kernels
  mask = 1.0_dp
  do iproc = myrank, (nproc-1), nrank

    print *, '#-- iproc=', iproc

    cijkl_kernel_sum = 0.0_dp
    do iker = 1, nkernel

      ! read kernel gll
      call sem_io_read_cijkl_kernel(kernel_dirs(iker), iproc, iregion, cijkl_kernel)

      ! read mask
      if (use_mask) then
        call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, mask_tag, mask)
      endif

      ! sum up kernels 
      cijkl_kernel_sum = cijkl_kernel_sum + cijkl_kernel*spread(mask, 1, 21)

    enddo ! iker

    ! write out cijkl kernel
    call sem_io_write_cijkl_kernel(out_dir, iproc, iregion, cijkl_kernel_sum)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
