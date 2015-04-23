! sem_inner_product <topo_dir> <kernel_dir1> <kernel_name1> 
!                   <kernel_dir2> <kernel_name2> 
!                   <kernel_out_dir> <kernel_out_name>
! 

subroutine selfdoc()

print *, 'Usage: sem_inner_product <topo_dir> <kernel_dir1> <kernel_name1> &
  & <kernel_dir2> <kernel_name2>' 
print *, 'out = inner product of (1, 2)'
print *, 'kernel_name: alpha_kernel, beta_kernel, vsv, mask_source etc.'
print *, 'kernel file: <kernel_dir>/proc000***_reg1_<kernel_name>.bin'

stop

end subroutine

program xsem_inner_product

  use sem_math

  implicit none

  !---- define variables
  integer :: ier, i
  integer, parameter :: nargs = 5
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

  read(args(1),'(a)') topo_dir 
  read(args(2),'(a)') kernel_dirs(1)
  read(args(3),'(a)') kernel_names(1)
  read(args(4),'(a)') kernel_dirs(2)
  read(args(5),'(a)') kernel_names(2)

  !---- program starts here
  L_COMPUTE_KERNEL_INTEGRAL = .true.
  L_WRITE_KERNEL_OUT = .false.
  opname = 'mul'
  call sem_math_op()

  print *, 'int(kernel)= ', kernel_out_integral,' volume= ',volume_integral
  print *, 'int(kernel)/volume= ',kernel_out_integral/volume_integral

end program xsem_inner_product
