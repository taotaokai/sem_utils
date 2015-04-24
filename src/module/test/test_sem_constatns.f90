program test

  use sem_constants

  implicit none

  integer :: iregion

  iregion = 1
  call sem_constants_set(iregion)

  print *, "NSPEC=", NSPEC, " NGLOB=", NGLOB

  print *, "NGLLX=", NGLLX

end program test
