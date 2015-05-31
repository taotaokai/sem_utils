module sem_tomography

!----------------------------------
! specification part
!----------------------------------

  use sem_constants
  use sem_io

  implicit none

private

  ! public subroutines
  !public :: sem_tomography_sum_kernels
  !public :: sem_tomography_sum_kernels_cijkl
  !public :: sem_tomography_reduce_kernel_cijkl_to_iso
  public :: sem_tomography_inner_product

!----------------------------------
! implementation part
!----------------------------------
contains

!!//////////////////////////
!! sum kernels that only consist of only one model property, e.g. rho_kernel
!subroutine sem_tomo_sum_kernels(kernel_sum, nker, kernel_dirs, iproc, iregion, kernel_name, flag_mask, mask_name)
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_sum 
!  integer, intent(in) :: nker, iproc, iregion
!  character(len=MAX_STRING_LEN), dimension(:), intent(in) :: kernel_dirs
!  character(len=*) :: kernel_name, mask_name
!  logical, intent(in) :: flag_mask
!
!  integer :: iker
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_mask
!  
!  kernel_sum = 0._CUSTOM_REAL
!
!  do iker = 1, nker
!
!    call sem_read_model(kernel, kernel_dirs(iker), iproc, iregion, kernel_name)
!
!    if (flag_mask) then
!      call sem_read_model(kernel_mask, kernel_dirs(iker), iproc, iregion, mask_name)
!    end if
!
!    kernel_sum = kernel_sum + kernel * kernel_mask
!
!  end do
!
!end subroutine
!
!!//////////////////////////
!! sum cijkl_kernels
!subroutine sem_sum_kernels_cijkl(kernel_cijkl_sum, nker, kernel_dirs, iproc, iregion, flag_mask, mask_name)
!
!  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_cijkl_sum
!  integer, intent(in) :: nker, iproc, iregion
!  character(len=MAX_STRING_LEN), dimension(:), intent(in) :: kernel_dirs
!  logical, intent(in) :: flag_mask
!  character(len=*) :: mask_name
!
!  integer :: iker, i
!  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_cijkl
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_mask
!  
!  kernel_cijkl_sum = 0._CUSTOM_REAL
!
!  do iker = 1, nker
!
!    call sem_read_kernel_cijkl(kernel_cijkl, kernel_dirs(iker), iproc, iregion)
!
!    if (flag_mask) then
!      call sem_read_model(kernel_mask, kernel_dirs(iker), iproc, iregion, mask_name)
!    end if
!
!    do i = 1, 21
!      kernel_cijkl_sum(i,:,:,:,:) = kernel_cijkl_sum(i,:,:,:,:) &
!                                   + kernel_cijkl(i,:,:,:,:) * kernel_mask
!    end do
!
!  end do
!
!end subroutine
!
!!//////////////////////////
!! reduce cijkl_kernel to lamda_2mu, mu kernels
!subroutine sem_reduce_kernel_cijkl_to_iso(kernel_lamda_2mu, kernel_mu, kernel_cijkl)
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_lamda_2mu
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: kernel_mu
!  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: kernel_cijkl
!
!  ! reduce cijkl_kernel to iso_kernel
!  kernel_lamda_2mu = kernel_cijkl(1,:,:,:,:) + &
!                     kernel_cijkl(2,:,:,:,:) + &
!                     kernel_cijkl(3,:,:,:,:) + &
!                     kernel_cijkl(7,:,:,:,:) + &
!                     kernel_cijkl(8,:,:,:,:) + &
!                     kernel_cijkl(12,:,:,:,:)
!
!  kernel_mu = -2._CUSTOM_REAL * (kernel_cijkl(2,:,:,:,:) + &
!                                 kernel_cijkl(3,:,:,:,:) + &
!                                 kernel_cijkl(8,:,:,:,:)) + &
!              kernel_cijkl(16,:,:,:,:) + &
!              kernel_cijkl(19,:,:,:,:) + &
!              kernel_cijkl(21,:,:,:,:)
!
!end subroutine

!!//////////////////////////
!! compute inner product of two iso models
!subroutine sem_tomography_inner_product_1(gll_model_1, gll_model_2, &
!  gll_volume, inner_product)
!!-calculate inner product between two sem gll models 
!! inner_product = volume integral of (model1 * model2)
!
!  real(dp), dimension(:,:,:,:), intent(in) :: gll_model_1, gll_model_2
!  real(dp), dimension(:,:,:,:), intent(in) :: gll_volume
!
!  real(dp), intent(out) :: inner_product
!
!  call kahan_sum4(gll_model_1 * gll_model_2 * gll_volume, inner_product)
!
!end subroutine

!//////////////////////////
subroutine kahan_sum4(model_array, model_sum)
!-Kahan summation for a rank-4 array
  
  real(dp), dimension(:,:,:,:), intent(in) :: model_array
  real(dp), intent(out) :: model_sum

  integer :: i,j,k,ispec
  integer :: dim1, dim2, dim3, dim4
  real(dp) :: t, y, c

  !-- get the size of the model array
  dim1 = size(model_array, dim=1)
  dim2 = size(model_array, dim=2)
  dim3 = size(model_array, dim=3)
  dim4 = size(model_array, dim=4)

  !-- kahan summation
  c = 0.0_dp
  model_sum = 0.0_dp

  do i4 = 1, dim4
    do i3 = 1, dim3
      do i2 = 1, dim2
        do i1 = 1, dim1
          y = model_array(i1, i2, i3, i4) - c
          t = model_sum + y
          c = (t - model_sum) - y
          model_sum = t
        enddo
      enddo
    enddo
  enddo

end subroutine

end module sem_tomography
