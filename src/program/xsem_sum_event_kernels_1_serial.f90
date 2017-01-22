subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_1"
  print '(a)', "    - sum up event kernels"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_1 "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <kernel_dir_list> <kernel_name> "
  print '(a)', "    <use_mask> <mask_tag> "
  print '(a)', "    <out_dir> <out_name>"
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
  print '(a)', "  (string) mask_tag:  tags of mask files"
  print '(a)', "  (string) out_dir:  directory to write out summed kernel files"
  print '(a)', "  (string) out_name:  <out_dir>/proc*_reg1_<out_name>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. only run in serial"
  print '(a)', "  2. mask files: <kernel_dir>/proc***_reg1_<mask_tag>.bin"
  print '(a)', ""

end subroutine


!//////////////////////////////////////////////////////////////////////////////
program xsem_sum_event_kernels_cijkl

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir_list
  character(len=MAX_STRING_LEN) :: kernel_name
  logical :: use_mask
  character(len=MAX_STRING_LEN) :: mask_tag 
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
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

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] check your inputs."
    stop
  endif

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
      print *, '[ERROR]: use_mask must be 0 or 1'
      stop
  end select
  read(args(6), '(a)') mask_tag
  read(args(7), '(a)') out_dir
  read(args(8), '(a)') out_name

  !====== read kernel_dir_list
  call sem_utils_read_line(kernel_dir_list, kernel_dirs, nkernel)

  !====== loop model slices 

  ! get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! initialize gll arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kernel_sum(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mask(NGLLX,NGLLY,NGLLZ,nspec))

  ! combine event kernels
  mask = 1.0_dp
  do iproc = 0, (nproc-1)

    print *, '====== iproc ', iproc

    kernel_sum = 0.0_dp
    do iker = 1, nkernel

      ! read kernel gll
      call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, kernel_name, kernel)

      ! read mask
      if (use_mask) then
        call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, mask_tag, mask)
      endif

      ! sum up kernels
      kernel_sum = kernel_sum + kernel * mask

    enddo ! iker

    ! write out gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, kernel_sum)

  enddo ! iproc

end program
