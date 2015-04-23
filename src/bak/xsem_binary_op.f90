! sem_binary_op <op_name> <model_dir_1>   <model_name_1> 
!                         <model_dir_2>   <model_name_2> 
!                         <model_dir_out> <model_name_out>

subroutine selfdoc()
  print *, 'Usage: sem_binary_op <opname> <model_dir_1> <model_name_1> &
    & <model_dir_2> <model_name_2> <model_dir_out> <model_name_out>' 
  print *, 'model_out = model_1 op model_2'
  print *, 'op_name: add, sub, mul, div'
  print *, 'model_name: vpv, rho_model, mu_model etc.'
  print *, 'model file read: <model_dir>/proc000***_reg1_<model_name>.bin'
  print *, 'model file write: <model_dir_out>/proc000***_reg1_<model_name_out>.bin'
  stop
end subroutine

program sem_binary_op

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use sem_IO

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 7

  !---- cmd line args
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: op_name
  character(len=MAX_STRING_LEN) :: model_dir_1, model_name_1 
  character(len=MAX_STRING_LEN) :: model_dir_2, model_name_2 
  character(len=MAX_STRING_LEN) :: model_dir_out, model_name_out

  !---- local variables
  integer :: iproc, i, ier
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: model_1, model_2, model_out 

  ! ============ program starts here =====================

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo

  read(args(1),'(a)') op_name
  read(args(2),'(a)') model_dir_1
  read(args(3),'(a)') model_name_1
  read(args(4),'(a)') model_dir_2
  read(args(5),'(a)') model_name_2
  read(args(6),'(a)') model_dir_out
  read(args(7),'(a)') model_name_out

  !---- initialize arrays
  call sem_set_dimension()
  allocate(  model_1(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_2(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_out(NGLLX,NGLLY,NGLLZ,NSPEC))
  
  !---- operate on each slice 
  do iproc = 0, NPROCTOT_VAL-1

    ! read models
    call sem_read_model(model_1, model_dir_1, iproc, iregion, model_name_1)
    call sem_read_model(model_2, model_dir_2, iproc, iregion, model_name_2)

    ! math operations
    select case (trim(op_name))
      case ('add')
        model_out = model_1 + model_2
      case ('sub')
        model_out = model_1 - model_2
      case ('mul')
        model_out = model_1 * model_2
      case ('div')
        model_out = model_1 / model_2
      ! default
      case default
        print *, 'Error invalid operation name: ', trim(opname)
        stop
    end select

    ! write out model
    call sem_write_model(model_out, model_dir_out, iproc, iregion, model_name_out)
  end do

end program
