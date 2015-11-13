subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_cijkl_to_lamda_mu "
  print '(a)', "    - reduce cijkl kernel to (lamda,mu) kernel"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_cijkl_to_lamda_mu \"
  print '(a)', "    <nproc> <kernel_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  "
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_cijkl_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for lamda,mu_kerenl"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_cijkl_to_lamda_mu

  use sem_constants
  use sem_io
  use sem_utils

  implicit none

  !===== declare variables
  ! command line args
  integer, parameter :: nargs = 3
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! kernel gll 
  real(dp), allocatable :: gll_cijkl(:,:,:,:,:)
  real(dp), allocatable :: gll_lamda(:,:,:,:), gll_mu(:,:,:,:)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_kernel_cijkl_to_lamda_mu: check your inputs."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), *) kernel_dir
  read(args(3), *) out_dir 

  !====== loop model slices 
  do iproc = 0, (nproc-1)

    ! read kernel gll
    call sem_io_read_gll_file_cijkl(kernel_dir, iproc, iregion, gll_cijkl)

    ! reduce cijkl kernel to lamda,mu kernel 
    kernel_lamda = kernel_cijkl(1,:,:,:,:) + &
                   kernel_cijkl(2,:,:,:,:) + &
                   kernel_cijkl(3,:,:,:,:) + &
                   kernel_cijkl(7,:,:,:,:) + &
                   kernel_cijkl(8,:,:,:,:) + &
                   kernel_cijkl(12,:,:,:,:)

    kernel_mu = kernel_cijkl(1,:,:,:,:)  + &
                kernel_cijkl(7,:,:,:,:)  + &
                kernel_cijkl(12,:,:,:,:) + &
                kernel_cijkl(16,:,:,:,:) + &
                kernel_cijkl(19,:,:,:,:) + &
                kernel_cijkl(21,:,:,:,:)

    ! write out lamda,mu kernel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'lamda_kerenl', gll_lamda)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'mu_kerenl', gll_mu)

  enddo ! iproc
