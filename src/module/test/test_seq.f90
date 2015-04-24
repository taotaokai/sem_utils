program test

  implicit none

  integer :: i, a(6)

  a = (/(i, i=1,6)/)

  print *, "a=", a

end program