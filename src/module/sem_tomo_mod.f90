module sem_tomo

!----------------------------------
! specification part
!----------------------------------

  use sem_constants
  use sem_io

  implicit none

private

  ! public subroutines
  public :: sem_tomo_sum_kernels
  public :: sem_tomo_sum_kernels_cijkl
  public :: sem_tomo_reduce_kernel_cijkl_to_iso
  public :: sem_tomo_inner_prod

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
! sum kernels that only consist of only one model property, e.g. rho_kernel
subroutine sem_tomo_sum_kernels(kernel_sum, nker, kernel_dirs, iproc, iregion, kernel_name, flag_mask, mask_name)

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
! compute inner product of two iso models
subroutine sem_inner_prod(mesh_data, nmodel, model_array1, model_array2, inner_prod)

  integer, intent(in) :: nmodel
  type(mesh), intent(in) :: mesh_data
  real(kind=CUSTOM_REAL), dimension(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: model_array1, model_array2
  real(kind=CUSTOM_REAL), intent(out) :: inner_prod

  integer :: i
  real(kind=8), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model_prod, volume_gll
  real(kind=8) :: sumval, model_sum

  call sem_mesh_get_volume_gll(mesh_data, volume_gll)

  model_sum = 0.d0
  do i = 1, nmodel
    model_prod = DBLE(model_array1(i,:,:,:,:)) &
                * DBLE(model_array2(i,:,:,:,:)) &
                * volume_gll
    call kahan_sum(model_prod, sumval)
    model_sum = model_sum + sumval 
  end do

  inner_prod = real(model_sum, kind=CUSTOM_REAL)

end subroutine

!//////////////////////////
subroutine kahan_sum(model_array, model_sum)
  
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

end module sem_tomo
