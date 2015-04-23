! sum_kernels_cijkl_to_iso <kernel_dir_list> <out_dir> <flag_mask> <mask_name>

subroutine selfdoc()
  print *, 'Usage: sum_kernels_cijkl_to_iso <kernel_dir_list> <out_dir> <flag_mask=0/1> <mask_name>'
  print *, 'files read: <kernel_dir>/proc000***_reg1_[cijkl,rho]_kernel.bin'
  print *, 'files written: <out_dir>/proc000***_reg1_[mu,lamda_2mu,rho]_kernel.bin'
  stop
end subroutine

program sum_kernels_cijkl_to_iso 

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use sem_IO
  use sem_tomography

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 4

  !---- local variables
  integer :: iproc, nker, i, ier
  logical :: flag_mask
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: kernel_dir_list, out_dir, model_name, mask_name
  character(len=MAX_STRING_LEN), allocatable :: kernel_dirs(:)
  real(CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: kernel_cijkl
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kernel_rho
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kernel_mu, kernel_lamda_2mu

  !---- get command line arguments 
  do i = 1, nargs 
    call get_command_argument(i,args(i), status=ier)
    if ((i <= 3) .and. (trim(args(i)) == '')) then
      call selfdoc()
    endif
  enddo
  read(args(1),'(a)') kernel_dir_list
  read(args(2),'(a)') out_dir
  select case (args(3))
    case ('0') 
      flag_mask = .false.
    case ('1') 
      flag_mask = .true.
    case default
      write(*,*) 'Error: <flag_mask> must be 0 or 1'
      call selfdoc() 
  end select
  if (flag_mask) then
    if (trim(args(4)) == '') then
      print *, 'Error: mask_name cannot be empty!'
      call selfdoc() 
    end if
    read(args(4),'(a)') mask_name
  end if

  !---- program starts here

  call read_text_list(kernel_dir_list, kernel_dirs, nker)
  call sem_set_dimension(iregion)

  allocate(kernel_cijkl(21,NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_rho(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_lamda_2mu(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_mu(NGLLX,NGLLY,NGLLZ,NSPEC))

  do iproc = 0, NPROCTOT_VAL-1

    ! sum kernels
    ! cijkl_kernel
    call sem_sum_kernels_cijkl(kernel_cijkl, nker, kernel_dirs, iproc, iregion, &
                               flag_mask, mask_name)
    call sem_reduce_kernel_cijkl_to_iso(kernel_lamda_2mu, kernel_mu, kernel_cijkl)

    model_name = 'rho_kernel'
    call sem_sum_kernels(kernel_rho, nker, kernel_dirs, iproc, iregion, & 
                         model_name, flag_mask, mask_name)

    ! write kernel
    model_name = 'mu_kernel'
    call sem_write_model(kernel_mu, out_dir, iproc, iregion, model_name)

    model_name = 'lamda_2mu_kernel'
    call sem_write_model(kernel_lamda_2mu, out_dir, iproc, iregion, model_name)

    model_name = 'rho_kernel'
    call sem_write_model(kernel_rho, out_dir, iproc, iregion, model_name)

  enddo

end program sum_kernels_cijkl_to_iso 
