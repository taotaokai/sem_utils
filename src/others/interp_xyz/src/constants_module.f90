module constants

  include "constants.h"

  include "values_from_mesher.h"

  real(CUSTOM_REAL), parameter :: XGLL(NGLLX) = & 
    (/-1.0,-sqrt(21.0)/7,0.0,sqrt(21.0)/7,1.0/)

  real(CUSTOM_REAL), parameter :: MESH_SIZE = &
    SNGL( (ANGULAR_WIDTH_XI_IN_DEGREES_VAL * DEGREES_TO_RADIANS/NEX_XI_VAL) * &
          R_UNIT_SPHERE )

end module constants
