subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_1"
  print '(a)', "    - sum up event kernels"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_sum_event_kernels_1 \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir_list> <kernel_name> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir_list:  list of event kernel directories (proc*_reg1_cijkl_kernel.bin)"
  print '(a)', "  (string) kernel_name:  kernel name to be summed (e.g. rho_kernel)"
  print '(a)', "  (string) out_dir:  output directory of summed kernel files"
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
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir_list
  character(len=MAX_STRING_LEN) :: kernel_name
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! list of kernel directories
  character(len=MAX_STRING_LEN), allocatable :: kernel_dirs(:)
  integer :: iker, nkernel

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec

  ! kernel gll 
  real(dp), allocatable :: gll(:,:,:,:)
  real(dp), allocatable :: gll_sum(:,:,:,:)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_sum_event_kernels_1: check your inputs."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), *) mesh_dir
  read(args(3), *) kernel_dir_list
  read(args(4), *) kernel_name
  read(args(5), *) out_dir 

  !====== read kernel_dir_list
  call sem_utils_read_line(kernel_dir_list, kernel_dirs, nkernel)

  !====== loop model slices 

  ! get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! initialize gll arrays 
  allocate(gll(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll_sum(NGLLX,NGLLY,NGLLZ,nspec))

  ! combine event kernels
  gll_sum = 0.0_dp
  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    do iker = 1, nkernel

      ! read kernel gll
      call sem_io_read_gll_file_1(kernel_dirs(iker), iproc, iregion, kernel_name, gll)

      gll_sum = gll_sum + gll

    enddo ! iker

    ! write out lamda,mu kernel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, kernel_name, gll_sum)

  enddo ! iproc

end program
