! xsem_gather_kernels: operates on several kernels and produce one kernel

subroutine selfdoc()

print *, 'Usage: sem_gather_kernels <opname> <kernel_dir_list> <kernel_name> &
  & <kernel_out_dir> <kernel_out_name>'
print *, 'opname: sum, prod'
print *, 'kernel_name: alpha_kernel, beta_kernel, vsv, mask_source etc.'
print *, 'kernel_dir_list: a list of kernel directories'
print *, 'kernel files: <kernel_dir>/proc000***_reg1_<kernel_name>.bin'

stop

end subroutine

program xsem_gather_kernels

  use sem_math

  implicit none

  !---- define variables
  integer :: ier, i
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs), dummy
  character(len=MAX_STRING_LEN) :: kernel_name, kernel_dir_list

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo

  read(args(1),'(a)') opname 
  read(args(2),'(a)') kernel_dir_list
  read(args(3),'(a)') kernel_name
  read(args(4),'(a)') kernel_out_dir
  read(args(5),'(a)') kernel_out_name

  !---- get kernel_dirs(:)
  open(unit=IIN, file=trim(kernel_dir_list), status='old', iostat=ier)
  if (ier /= 0) then
     print *, 'Error open file: ', trim(kernel_dir_list)
     stop
  endif
  kernel_num = 0 ! get number of kernel_dirs 
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    kernel_num = kernel_num+1
  enddo
  allocate(kernel_dirs(kernel_num),kernel_names(kernel_num))
  rewind(IIN)
  do i = 1, kernel_num
    read(IIN,'(a)',iostat=ier) kernel_dirs(i)
  enddo
  close(IIN)
  kernel_names = kernel_name

  !---- program starts here
  L_COMPUTE_KERNEL_INTEGRAL = .false.
  L_WRITE_KERNEL_OUT = .true.
  call sem_math_op()

end program xsem_gather_kernels
