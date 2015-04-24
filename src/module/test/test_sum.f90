program test

  implicit none

  integer :: a(6), b(6)

  a = [1, 2, 3, 4, 5, 6]
  b = [1, 2, 3, 4, 5, 6]

  print *, "sum(a+b) = ", sum(a+b, DIM=2)

end program