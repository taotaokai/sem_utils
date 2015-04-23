! report kernel info: max,min,mean

subroutine selfdoc()

print *, 'Usage: sem_check <kernel_dir> <kernel_name>' 
print *, 'kernel_name: alpha_kernel, beta_kernel, vsv, mask_source etc.'
print *, 'kernel file: <kernel_dir>/proc000***_reg1_<kernel_name>.bin'

stop

end subroutine

program xsem_check

  use sem_math

  implicit none

  !---- define variables
  integer :: ier, i
  integer, parameter :: nargs = 2
  character(len=MAX_STRING_LEN) :: args(nargs)

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo

  kernel_num = 1
  allocate(kernel_dirs(kernel_num), kernel_names(kernel_num))

  read(args(1),'(a)') kernel_dirs(1)
  read(args(2),'(a)') kernel_names(1)

  !---- program starts here
  call sem_check()

end program xsem_check
