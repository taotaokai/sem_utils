module sem_data

!----------------------------------
! specification part
!----------------------------------

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE

  implicit none

  ! inputs
  integer :: myrank, iregion
  character(len=MAX_STRING_LEN) :: MESH_DIR, MODEL_DIR, KERNEL_DIR

  ! dimension
  integer :: NSPEC, NGLOB

  ! mesh
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x, y, z
  integer, dimension(:,:,:,:),allocatable :: ibool
  logical, dimension(:),allocatable :: ispec_is_tiso

  ! jacobian of element shape function

  ! kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: kernel_cijkl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_rho

  ! kernel_mask
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_mask

  ! model
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    model_vpv, model_vph, model_vsv, model_vsh, model_eta, model_rho 

private
 
  character(len=MAX_STRING_LEN) :: basename
  logical :: array_dimension_set = .false.
  logical :: mesh_array_allocated = .false.
  logical :: model_array_allocated = .false.
  logical :: kernel_array_allocated = .false.
  logical :: mask_array_allocated = .false.

!---- public routines
  public :: read_model, read_kernel_cijkl_rho, read_kernel_mask

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
subroutine set_dimension()
  select case (iregion)
    case (1)
      NSPEC = NSPEC_CRUST_MANTLE
      NGLOB = NGLOB_CRUST_MANTLE
!   case (2)
!     NSPEC = NSPEC_CRUST_MANTLE
!     NGLOB = NGLOB_CRUST_MANTLE
!   case (3)
!     NSPEC = NSPEC_CRUST_MANTLE
!     NGLOB = NGLOB_CRUST_MANTLE
    case default
      write(*,*) 'Error: unrecognized region code ( ', iregion, ' )'
      stop
  end select
  array_dimension_set = .true.
end subroutine

!//////////////////////////
subroutine cread_basename(basedir)

  character(len=MAX_STRING_LEN) :: basedir
  write(basename,"(a,'/proc',i6.6,'_reg',i1,'_')") trim(basedir),myrank,iregion

end subroutine

!//////////////////////////
subroutine read_mesh()

  call set_dimension()

  ! initialize arrays
  if (.not. mesh_array_allocated) then
    allocate(x(NGLOB),y(NGLOB),z(NGLOB),ispec_is_tiso(NGLOB))
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC))
  end if

  call cread_basename(MESH_DIR)

end subroutine

!//////////////////////////
subroutine read_model()

  call set_dimension()

  ! initialize model arrays

end subroutine

!//////////////////////////
subroutine read_kernel_cijkl_rho()

  call set_dimension()

  ! initialize kernel arrays

end subroutine

!//////////////////////////
subroutine read_kernel_mask()

  call set_dimension()

  ! initialize kernel arrays

end subroutine

end module sem_data
