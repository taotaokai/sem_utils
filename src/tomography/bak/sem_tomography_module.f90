module sem_tomography

!----------------------------------
! specification part
!----------------------------------

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA

  use sem_IO

  implicit none

private

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll
  ! array with all the integration weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  logical :: flag_gll_weight = .false.

  ! public subroutines
  public :: read_text_list
  public :: sem_sum_kernels
  public :: sem_sum_kernels_cijkl
  public :: sem_reduce_kernel_cijkl_to_iso
  public :: sem_volume_weight_gll
  public :: sem_inner_product_iso
  public :: timestamp

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
! read text file line by line
subroutine read_text_list(filename, slist, nline)

  character(len=MAX_STRING_LEN), intent(in) :: filename

  integer, intent(out) :: nline 
  character(len=MAX_STRING_LEN), dimension(:), allocatable, intent(out) :: slist

  integer :: i, ier
  character(len=MAX_STRING_LEN) :: dummy

  open(unit=IIN, file=trim(filename), status='old', iostat=ier)
  if (ier /= 0) then
     write(*,*) 'Error open file: ', trim(filename)
     stop
  endif
  nline = 0 ! get number of lines
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    nline = nline + 1
  enddo
  if (.not. allocated(slist)) then
    allocate(slist(nline))
  elseif (size(slist) /= nline) then
    deallocate(slist)
    allocate(slist(nline))
  end if
  rewind(IIN)
  do i = 1, nline 
    read(IIN,'(a)',iostat=ier) slist(i)
  enddo
  close(IIN)

end subroutine

!//////////////////////////
! sum kernels that only consist of only one model property, e.g. rho_kernel
subroutine sem_sum_kernels(kernel_sum, nker, kernel_dirs, iproc, iregion, kernel_name, flag_mask, mask_name)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_sum 
  integer, intent(in) :: nker, iproc, iregion
  character(len=MAX_STRING_LEN), dimension(:), intent(in) :: kernel_dirs
  character(len=*) :: kernel_name, mask_name
  logical, intent(in) :: flag_mask

  integer :: iker
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_mask
  
  kernel_sum = 0._CUSTOM_REAL

  do iker = 1, nker

    call sem_read_model(kernel, kernel_dirs(iker), iproc, iregion, kernel_name)

    if (flag_mask) then
      call sem_read_model(kernel_mask, kernel_dirs(iker), iproc, iregion, mask_name)
    end if

    kernel_sum = kernel_sum + kernel * kernel_mask

  end do

end subroutine

!//////////////////////////
! sum cijkl_kernels
subroutine sem_sum_kernels_cijkl(kernel_cijkl_sum, nker, kernel_dirs, iproc, iregion, flag_mask, mask_name)

  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_cijkl_sum
  integer, intent(in) :: nker, iproc, iregion
  character(len=MAX_STRING_LEN), dimension(:), intent(in) :: kernel_dirs
  logical, intent(in) :: flag_mask
  character(len=*) :: mask_name

  integer :: iker, i
  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_cijkl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_mask
  
  kernel_cijkl_sum = 0._CUSTOM_REAL

  do iker = 1, nker

    call sem_read_kernel_cijkl(kernel_cijkl, kernel_dirs(iker), iproc, iregion)

    if (flag_mask) then
      call sem_read_model(kernel_mask, kernel_dirs(iker), iproc, iregion, mask_name)
    end if

    do i = 1, 21
      kernel_cijkl_sum(i,:,:,:,:) = kernel_cijkl_sum(i,:,:,:,:) &
                                   + kernel_cijkl(i,:,:,:,:) * kernel_mask
    end do

  end do

end subroutine

!//////////////////////////
! reduce cijkl_kernel to lamda_2mu, mu kernels
subroutine sem_reduce_kernel_cijkl_to_iso(kernel_lamda_2mu, kernel_mu, kernel_cijkl)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_lamda_2mu
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_mu
  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: kernel_cijkl

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

end subroutine

!//////////////////////////
! compute the integration weight for gll quadrature in cube
subroutine compute_gll_weight()

  integer :: i,j,k

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

  flag_gll_weight = .true.

end subroutine

!//////////////////////////
! compute volume integration weight on gll points in each spectral element
! wgll_vol = jacobian * wgll_cube
! jacobian = det(d(x,y,z)/d(xi,eta,gamma)), distortion from cube to actual element
subroutine sem_volume_weight_gll(mesh_data, wgll_vol)

  type (mesh), intent(in) :: mesh_data
  real(kind=8), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: wgll_vol 

  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  ! calculate gll quad weight in the reference cube
  if (.not. flag_gll_weight) then
    call compute_gll_weight()
  end if

  ! calculate volume integration weight on gll
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! gets derivatives of ux, uy and uz with respect to x, y and z
             xixl = mesh_data%xix(i,j,k,ispec)
             xiyl = mesh_data%xiy(i,j,k,ispec)
             xizl = mesh_data%xiz(i,j,k,ispec)
            etaxl = mesh_data%etax(i,j,k,ispec)
            etayl = mesh_data%etay(i,j,k,ispec)
            etazl = mesh_data%etaz(i,j,k,ispec)
          gammaxl = mesh_data%gammax(i,j,k,ispec)
          gammayl = mesh_data%gammay(i,j,k,ispec)
          gammazl = mesh_data%gammaz(i,j,k,ispec)

          ! computes the integration weight 
          wgll_vol(i,j,k,ispec) = wgll_cube(i,j,k) / &
                       (  xixl*(etayl*gammazl-etazl*gammayl)  &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl)  &
                        + xizl*(etaxl*gammayl-etayl*gammaxl) )
        enddo
      enddo
    enddo
  enddo

end subroutine

!//////////////////////////
! compute inner product between two iso models (lamda_2mu, mu, rho)
subroutine sem_inner_product_iso(model_iso1, model_iso2, wgll_vol, inner_product)

  type(model_iso), intent(in) :: model_iso1, model_iso2
  real(kind=8), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: wgll_vol
  real(kind=8), intent(out) :: inner_product 

  real(kind=8), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model_prod

  model_prod = wgll_vol * ( &
    model_iso1%lamda_2mu * model_iso2%lamda_2mu + &
    model_iso1%mu * model_iso2%mu + &
    model_iso1%rho * model_iso2%rho )

  call Kahan_sum(model_prod, inner_product)

end subroutine

!//////////////////////////
subroutine Kahan_sum(model_array, model_sum)
  
  real(kind=8), dimension(:,:,:,:), intent(in) :: model_array
  real(kind=8), intent(out) :: model_sum

  integer :: i,j,k,ispec
  real(kind=8) :: t, y, c

  c = 0.d0
  model_sum = 0.d0

  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          y = model_array(i,j,k,ispec) - c
          t = model_sum + y
          c = (t - model_sum) - y
          model_sum = t
        end do
      end do
    end do
  end do

end subroutine

!//////////////////////////
! time stamp: yyyy-mm-ddThh:mm:ss.sss
character(len=MAX_STRING_LEN) function timestamp()

  implicit none
  
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  integer :: values(8)
  integer y, m, d, h, n, s, ms

  call date_and_time(date, time, zone, values)
  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  ms = values(8)

  write(timestamp,'(i4,"-",i2,"-",i2,"T",i2,":",i2,":",i2,".",i3)') &
    y, m, d, h, n, s, ms

  return

end function 


end module sem_tomography
