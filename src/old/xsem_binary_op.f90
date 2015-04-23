! sem_binary_op <opname> <kernel_dir1> <kernel_name1> 
!                        <kernel_dir2> <kernel_name2> 
!                        <kernel_out_dir> <kernel_out_name>
! 

subroutine selfdoc()

print *, 'Usage: sem_binary_op <opname> <kernel_dir1> <kernel_name1> &
  & <kernel_dir2> <kernel_name2> <kernel_out_dir> <kernel_out_name>' 
print *, 'out = 1 op 2'
print *, 'opname: add, sub, mul, div' 
print *, 'kernel_name: alpha_kernel, beta_kernel, vsv, mask_source etc.'
print *, 'kernel file: <kernel_dir>/proc000***_reg1_<kernel_name>.bin'

stop

end subroutine

program xsem_binary_op

  use sem_math

  implicit none

  !---- define variables
  integer :: ier, i
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo

  kernel_num = 2
  allocate(kernel_dirs(kernel_num), kernel_names(kernel_num))

  read(args(1),'(a)') opname
  read(args(2),'(a)') kernel_dirs(1)
  read(args(3),'(a)') kernel_names(1)
  read(args(4),'(a)') kernel_dirs(2)
  read(args(5),'(a)') kernel_names(2)
  read(args(6),'(a)') kernel_out_dir
  read(args(7),'(a)') kernel_out_name

  !---- program starts here
  L_COMPUTE_KERNEL_INTEGRAL = .false.
  L_WRITE_KERNEL_OUT = .true.
  call sem_math_op()

end program xsem_binary_op
