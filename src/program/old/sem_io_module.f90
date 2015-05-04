! math operations on kernel files

module sem_io

! ================ Specification Part ===============

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE

  implicit none

private

  !---- local variables
  integer :: myrank, ier, iker, i, j, k, ispec
  character(len=MAX_STRING_LEN) :: file_name

  !---- public arrays
  public :: MAX_STRING_LEN, IIN

  !---- public routines
  public :: sem_read_model, sem_write_model, sem_read_mesh

! ================ Implementation Part ==============
contains

subroutine compute_gll_weight()

  ! GLL points
  wgll_cube = 0.0d0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
      enddo
    enddo
  enddo

  L_GLL_WEIGHT_COMPUTED = .true.

end subroutine compute_gll_weight

!/////////////////////////////////////////////////!

subroutine read_kernels_myrank()
! read all kernels of {myrank} into one array

  implicit none

  ! initializes arrays
  if (.not. allocated(kernels)) then
    allocate(kernels(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,kernel_num), stat=ier)
    if (ier /= 0) stop 'Error allocating arrays (read_kernels_myrank)'
  endif

  ! read all kernels of myrank
  do iker = 1, kernel_num

    write(file_name, '(a,i6.6,a)') trim(kernel_dirs(iker))//'/proc', &
      myrank, trim(REG)//trim(kernel_names(iker))//'.bin'

    open(IIN,file=trim(file_name),status='old', &
      form='unformatted',action='read',iostat=ier)

    if (ier /= 0) then
      write(*,*) 'Error open file: ',trim(file_name)
      stop
    endif

    read(IIN) kernels(:,:,:,:,iker)

    close(IIN)
  enddo

end subroutine read_kernels_myrank

!/////////////////////////////////////////////////!

subroutine write_kernel_out_myrank()
! write kernel_out of myrank

  ! initializes arrays
  if (.not. allocated(kernel_out)) then
    stop 'Error: kernel_out not allocated (write_kernel_out_myrank)'
  endif

  write(file_name,'(a,i6.6,a)') trim(kernel_out_dir)//'/proc', &
    myrank, trim(REG)//trim(kernel_out_name)//'.bin'

  open(IOUT,file=trim(file_name),form='unformatted', &
    status='unknown',action='write',iostat=ier)

  if (ier /= 0) then
    write(*,*) 'Error open file: ',trim(file_name)
    stop
  endif

  write(IOUT) kernel_out

  close(IOUT)
end subroutine write_kernel_out_myrank

!/////////////////////////////////////////////////!

subroutine sem_math_op()
! math operations

  implicit none

  ! local variables
  real(CUSTOM_REAL) :: kl_integral_myrank, volume_myrank

  ! initialize
  if (L_COMPUTE_KERNEL_INTEGRAL) then
    kernel_out_integral = 0.0_CUSTOM_REAL
    volume_integral = 0.0_CUSTOM_REAL
  endif

  if (.not. allocated(kernel_out)) then
    allocate(kernel_out(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), stat=ier)
    if (ier /= 0) stop 'Error allocating arrays (sem_math_op:kernel_out)'
  endif

  ! loop each slice 
  do myrank = 0, NPROCTOT_VAL-1

    call read_kernels_myrank()

    ! math operations
    select case (trim(opname))
      ! binary operation: 2 to 1
      case ('add')
        kernel_out = kernels(:,:,:,:,1) + kernels(:,:,:,:,2)
      case ('sub')
        kernel_out = kernels(:,:,:,:,1) - kernels(:,:,:,:,2)
      case ('mul')
        kernel_out = kernels(:,:,:,:,1) * kernels(:,:,:,:,2)
      case ('div')
        kernel_out = kernels(:,:,:,:,1) / kernels(:,:,:,:,2)
      ! many to 1
      case ('sum')
        kernel_out = sum(kernels,dim=5)
      case ('prod')
        kernel_out = product(kernels,dim=5)
      ! default
      case default
        print *, 'Error invalid operation: ', trim(opname)
        stop
    endselect

    if (L_WRITE_KERNEL_OUT) then
      call write_kernel_out_myrank()
    endif

    if (L_COMPUTE_KERNEL_INTEGRAL) then
      call compute_kernel_integral(kernel_out, kl_integral_myrank, volume_myrank)
      kernel_out_integral = kernel_out_integral + kl_integral_myrank
      volume_integral = volume_integral + volume_myrank
      print *, 'myrank= ',myrank,' int(kernel)= ',kl_integral_myrank,' volume= ',volume_myrank
    endif

  enddo

end subroutine sem_math_op

!///////////////////////////////////////////!

! volume integral of kernel
subroutine compute_kernel_integral(kernel, kernel_integral, volume_integral)

  implicit none

  ! input
  real(CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
    intent(in) ::  kernel
  real(CUSTOM_REAL), intent(out) :: kernel_integral, volume_integral

  ! local vars
  real(CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: jacobian
  real(CUSTOM_REAL) :: volumel

  ! initialize
  if (.not. L_GLL_WEIGHT_COMPUTED) then
    call compute_gll_weight()
  endif

  ! builds jacobian
  call compute_jacobian(jacobian)

  ! compute volume integral
  kernel_integral = 0.0_CUSTOM_REAL
  volume_integral = 0.0_CUSTOM_REAL
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          volumel = jacobian(i,j,k,ispec)*wgll_cube(i,j,k)
          kernel_integral = kernel_integral + volumel*kernel(i,j,k,ispec)
          volume_integral = volume_integral + volumel
        enddo
      enddo
    enddo
  enddo

end subroutine compute_kernel_integral

!///////////////////////////////////////////!

subroutine compute_jacobian(jacobian)
! computes volume element associated with points

  implicit none

  !input
  real(CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
    intent(out) :: jacobian

  ! dummy
  integer :: ival
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: dummy
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: dummy_ibool
  integer, dimension(NSPEC_CRUST_MANTLE) :: dummy_idoubling
  logical, dimension(NSPEC_CRUST_MANTLE) :: dummy_ispec_is_tiso

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL) :: &
    xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, jacobianl
  integer :: i,j,k,ispec

  ! initializes
  jacobian(:,:,:,:) = 0.0_CUSTOM_REAL

  ! builds jacobian
  write(file_name,'(a,i6.6,a)') trim(topo_dir)//'/proc', myrank, trim(REG)//'solver_data.bin'
  open(IIN,file=trim(file_name),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error open file: ',trim(file_name)
    stop
  endif

  read(IIN) ival ! nspec
  read(IIN) ival ! nglob

  read(IIN) dummy ! x
  read(IIN) dummy ! y
  read(IIN) dummy ! z

  read(IIN) dummy_ibool
  read(IIN) dummy_idoubling
  read(IIN) dummy_ispec_is_tiso

  read(IIN) xix
  read(IIN) xiy
  read(IIN) xiz
  read(IIN) etax
  read(IIN) etay
  read(IIN) etaz
  read(IIN) gammax
  read(IIN) gammay
  read(IIN) gammaz
  close(IIN)

  ! calculates jacobian
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! gets derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          ! computes the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))
          jacobian(i,j,k,ispec) = jacobianl
        enddo
      enddo
    enddo
  enddo

end subroutine compute_jacobian

!///////////////////////////////////////////!

subroutine sem_check()

  implicit none

  ! loop each rank
  do myrank = 0, NPROCTOT_VAL-1

    print *, 'myrank=', myrank
    call read_kernels_myrank()

    ! loop each kernel files
    do iker = 1, kernel_num
      ! write out kernel info
      print *, trim(kernel_dirs(iker)),'  ',trim(kernel_names(iker))
      write(*,'(a,3E15.7)') 'min/max/mean', &
        minval(kernels(:,:,:,:,iker)), &
        maxval(kernels(:,:,:,:,iker)), &
        sum(kernels(:,:,:,:,iker))/size(kernels(:,:,:,:,iker))
    enddo

  enddo

end subroutine sem_check

end module sem_math
