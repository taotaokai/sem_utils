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
  ! define mesh
  public :: NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB
  ! define gll function
  public :: GAUSSALPHA, GAUSSBETA

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
      write(*,*) 'ERROR:sem_constants: unrecognized region code ( ', iregion, ' )'
      stop
  end select

end subroutine


end module sem_constants
