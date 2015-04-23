! sum_kernels_cijkl_to_iso <kernel_dir_list> <out_dir> <mask_source_flag>

subroutine selfdoc()
  print *, 'Usage: sum_kernels_cijkl_to_iso <kernel_dir_list> <out_dir> <mask_source_flag=0/1>'
  print *, 'files read: <kernel_dir>/proc000***_reg1_[cijkl,rho]_kernel.bin'
  print *, 'files written: <out_dir>/proc000***_reg1_[mu,lamda_2mu,rho]_kernel.bin'
  stop
end subroutine

program sum_kernels_cijkl_to_iso 

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,NPROCTOT_VAL

  implicit none

  !---- define variables
  integer :: ier, i
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 3
  character(len=MAX_STRING_LEN) :: args(nargs), dummy
  character(len=MAX_STRING_LEN) :: kernel_dir_list

  !---- local variables
  integer :: myrank, iker, NKERNEL 
  logical :: mask_source_flag
  character(len=MAX_STRING_LEN) :: file_name, out_dir
  character(len=MAX_STRING_LEN), allocatable :: kernel_dirs(:) 
  real(CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: kernels_in_cijkl
  real(CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    kernels_in_rho, kernels_in_mask
  real(CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    kernel_mu, kernel_lamda_2mu, kernel_rho

  !---- get command line arguments 
  do i = 1, nargs 
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo
  read(args(1),'(a)') kernel_dir_list
  read(args(2),'(a)') out_dir
  select case (args(3))
    case ('0') 
      mask_source_flag = .false.
    case ('1') 
      mask_source_flag = .true.
    case default 
      write(*,*) '<mask_source_flag> must be 0 or 1'
      stop
  end select

  !---- program starts here
  call read_list()
  do myrank = 0, NPROCTOT_VAL-1
    call read_model()
    call sum_kernels_iso()
    call write_model('mu_kernel',kernel_mu)
    call write_model('lamda_2mu_kernel',kernel_lamda_2mu)
    call write_model('rho_kernel',kernel_rho)
  enddo

!==============================
contains

subroutine read_list()

  open(unit=IIN, file=trim(kernel_dir_list), status='old', iostat=ier)
  if (ier /= 0) then
     write(*,*) 'Error open file: ', trim(kernel_dir_list)
     stop
  endif
  NKERNEL = 0 ! get number of kernel_dirs 
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    NKERNEL = NKERNEL + 1
  enddo
  allocate(kernel_dirs(NKERNEL))
  rewind(IIN)
  do i = 1, NKERNEL
    read(IIN,'(a)',iostat=ier) kernel_dirs(i)
  enddo
  close(IIN)

end subroutine

!//////////////////////////////!
subroutine read_model()

  ! initialize arrays
  if (.not. allocated(kernels_in_cijkl)) then
    allocate(kernels_in_cijkl(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,NKERNEL))
    allocate(kernels_in_mask(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,NKERNEL))
    if (mask_source_flag) then
      allocate(kernels_in_mask(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,NKERNEL))
    endif
  endif

  ! read kernels
  do iker = 1, NKERNEL

    ! cijkl_kernel
    write(file_name,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
      trim(kernel_dirs(i)),myrank,iregion,'cijkl_kernel'
    open(IIN,file=trim(file_name),status='old', &
      form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) 'Error open file for read: ', trim(file_name)
      stop
    endif
    read(IIN) kernels_in_cijkl(:,:,:,:,:,iker)
    close(IIN)

    ! rho_kernel
    write(file_name,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
      trim(kernel_dirs(i)),myrank,iregion,'rho_kernel'
    open(IIN,file=trim(file_name),status='old', &
      form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) 'Error open file for read: ', trim(file_name)
      stop
    endif
    read(IIN) kernels_in_rho(:,:,:,:,iker)
    close(IIN)

    ! mask_source
    if (mask_source_flag) then
      write(file_name,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
        trim(kernel_dirs(i)),myrank,iregion,'mask_source'
      open(IIN,file=trim(file_name),status='old', &
        form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) 'Error open file for read: ', trim(file_name)
        stop
      endif
      read(IIN) kernels_in_mask(:,:,:,:,iker)
      close(IIN)
    endif
  enddo

end subroutine

!//////////////////////////////!
! reduce cijkl_kernel to iso_kernel
! now only plain summation with source mask
subroutine sum_kernels_iso()

  if (mask_source_flag) then

    kernel_lamda_2mu = sum( &
      (kernels_in_cijkl(1,:,:,:,:,:) + &
       kernels_in_cijkl(2,:,:,:,:,:) + &
       kernels_in_cijkl(3,:,:,:,:,:) + &
       kernels_in_cijkl(7,:,:,:,:,:) + &
       kernels_in_cijkl(8,:,:,:,:,:) + &
       kernels_in_cijkl(12,:,:,:,:,:)) * kernels_in_mask, dim=5)

    kernel_mu = sum( &
      (-2._CUSTOM_REAL * (kernels_in_cijkl(2,:,:,:,:,:) + &
                          kernels_in_cijkl(3,:,:,:,:,:) + &
                          kernels_in_cijkl(8,:,:,:,:,:)) + &
       kernels_in_cijkl(16,:,:,:,:,:) + &
       kernels_in_cijkl(19,:,:,:,:,:) + &
       kernels_in_cijkl(21,:,:,:,:,:)) * kernels_in_mask, dim=5)

    kernel_rho = sum(kernels_in_rho * kernels_in_mask, dim=5) 

  else

    kernel_lamda_2mu = sum( &
      (kernels_in_cijkl(1,:,:,:,:,:) + &
       kernels_in_cijkl(2,:,:,:,:,:) + &
       kernels_in_cijkl(3,:,:,:,:,:) + &
       kernels_in_cijkl(7,:,:,:,:,:) + &
       kernels_in_cijkl(8,:,:,:,:,:) + &
       kernels_in_cijkl(12,:,:,:,:,:)), dim=5)

    kernel_mu = sum( &
      (-2._CUSTOM_REAL * (kernels_in_cijkl(2,:,:,:,:,:) + &
                          kernels_in_cijkl(3,:,:,:,:,:) + &
                          kernels_in_cijkl(8,:,:,:,:,:)) + &
       kernels_in_cijkl(16,:,:,:,:,:) + &
       kernels_in_cijkl(19,:,:,:,:,:) + &
       kernels_in_cijkl(21,:,:,:,:,:)), dim=5)

    kernel_rho = sum(kernels_in_rho, dim=5) 

  endif

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

end program sum_kernels_cijkl_to_iso 
