! kernel_cijkl_to_iso <in_dir> <out_dir>
! 

subroutine selfdoc()
  print *, 'Usage: sem_kernel_cijkl_to_iso <in_dir> <out_dir>'
  print *, 'files read: <in_dir>/proc000***_reg1_cijkl_kernel.bin'
  print *, 'files written: <out_dir>/proc000***_reg1_<mu/lamda_2mu>_kernel.bin'
  stop
end subroutine

program kernel_cijkl_to_iso

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,NPROCTOT_VAL

  implicit none

  !---- define variables
  integer :: ier, i
  integer, parameter :: nargs = 2
  integer, parameter :: iregion = 1
  character(len=MAX_STRING_LEN) :: args(nargs)

  !---- local variables
  integer :: myrank
  character(len=MAX_STRING_LEN) :: file_name, in_dir, out_dir
  real(CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    kernel_cijkl
  real(CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    kernel_mu, kernel_lamda_2mu

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo

  read(args(1),'(a)') in_dir
  read(args(2),'(a)') out_dir

  !---- program starts here
  do myrank = 0, NPROCTOT_VAL-1
    call read_model()

    ! reduce cijkl_kernel to iso_kernel
    kernel_lamda_2mu = kernel_cijkl(1,:,:,:,:) + &
                       kernel_cijkl(2,:,:,:,:) + &
                       kernel_cijkl(3,:,:,:,:) + &
                       kernel_cijkl(7,:,:,:,:) + &
                       kernel_cijkl(8,:,:,:,:) + &
                       kernel_cijkl(12,:,:,:,:)
    kernel_mu = -2._CUSTOM_REAL * (kernel_cijkl(2,:,:,:,:) + &
                                  kernel_cijkl(3,:,:,:,:) + &
                                  kernel_cijkl(8,:,:,:,:)) + &
                kernel_cijkl(16,:,:,:,:) + &
                kernel_cijkl(19,:,:,:,:) + &
                kernel_cijkl(21,:,:,:,:)

    call write_model('mu_kernel',kernel_mu)
    call write_model('lamda_2mu_kernel',kernel_lamda_2mu)
  enddo

!==============================
contains

subroutine read_model()

  write(file_name,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(in_dir),myrank,iregion,'cijkl_kernel'
  open(IIN,file=trim(file_name),status='old', &
    form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error open file for read: ', trim(file_name)
    stop
  endif
  read(IIN) kernel_cijkl(:,:,:,:,:)
  close(IIN)

end subroutine

!//////////////////////////////!
subroutine write_model(kernel_name, kernel_array)

  character(len=*), intent(in) :: kernel_name 
  real(CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: kernel_array

  write(file_name,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(out_dir),myrank,iregion,trim(kernel_name)
  open(IOUT,file=trim(file_name),form='unformatted', &
    status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error open file for write: ', trim(file_name)
    stop
  endif
  write(IOUT) kernel_array
  close(IOUT)

end subroutine

end program kernel_cijkl_to_iso
