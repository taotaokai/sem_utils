module sem_constants

!----------------------------------
! specification part
!----------------------------------

  implicit none

  private

  include "constants.h"

  include "values_from_mesher.h"

! include 'constants_tomography.h'

  ! array dimension
  integer :: NSPEC, NGLOB

  !---- public varaibles
  public :: CUSTOM_REAL, MAX_STRING_LEN, IIN, IOUT
  public :: NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB
  ! sem_mesh: 
  public :: GAUSSALPHA, GAUSSBETA
  public :: NDIM, NGNOD
  public :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL
  public :: DEGREES_TO_RADIANS, R_UNIT_SPHERE
  public :: NEX_XI_VAL
  public :: NUM_ITER

  !---- public operations
  public :: sem_constants_set

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
subroutine sem_constants_set(iregion)

  integer, intent(in) :: iregion

  select case (iregion)
    case (IREGION_CRUST_MANTLE)
      NSPEC = NSPEC_CRUST_MANTLE
      NGLOB = NGLOB_CRUST_MANTLE
!   case (IREGION_OUTER_CORE)
!     NSPEC = NSPEC_CRUST_MANTLE
!     NGLOB = NGLOB_CRUST_MANTLE
!   case (IREGION_INNER_CORE)
!     NSPEC = NSPEC_CRUST_MANTLE
!     NGLOB = NGLOB_CRUST_MANTLE
    case default
      write(*,*) 'sem_constants_set: unrecognized region code ( ', iregion, ' )'
      stop
  end select

end subroutine


end module sem_constants
