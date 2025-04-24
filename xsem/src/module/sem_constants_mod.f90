module sem_constants

!----------------------------------
! specification part
!----------------------------------

  implicit none

  include "constants.h"

  include "values_from_mesher.h"

  ! real kind
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: sp = kind(1.0)

! include 'constants_tomography.h'

  ! array dimension
  !integer :: NSPEC, NGLOB

  !---- public varaibles
  !public :: CUSTOM_REAL, MAX_STRING_LEN, IIN, IOUT
  !public :: NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB
  ! sem_mesh: 
  !public :: GAUSSALPHA, GAUSSBETA
  !public :: NDIM, NGNOD
  !public :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL
  !public :: DEGREES_TO_RADIANS, R_UNIT_SPHERE
  !public :: NEX_XI_VAL
  !public :: NUM_ITER

  !---- public operations
  !public :: sem_constants_set

!----------------------------------
! implementation part
!----------------------------------
!contains

!//////////////////////////
!subroutine sem_constants_getdim(iregion, nspec, nglob)
!
!  integer, intent(in) :: iregion
!
!  integer, intent(out) :: iregion
!
!  select case (iregion)
!    case (IREGION_CRUST_MANTLE)
!      nspec = NSPEC_CRUST_MANTLE
!      nglob = NGLOB_CRUST_MANTLE
!!   case (IREGION_OUTER_CORE)
!!     nspec = NSPEC_CRUST_MANTLE
!!     nglob = NGLOB_CRUST_MANTLE
!!   case (IREGION_INNER_CORE)
!!     nspec = NSPEC_CRUST_MANTLE
!!     nglob = NGLOB_CRUST_MANTLE
!    case default
!      write(*,*) '[ERROR] sem_constants_getdim: unrecognized region code ( ', iregion, ' )'
!      stop
!  end select
!
!end subroutine


end module sem_constants
