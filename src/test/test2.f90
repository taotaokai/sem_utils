module test2

  use test

contains

subroutine tryit()

  integer, dimension(n) :: x

  print *, 'size(x)=', size(x)

end subroutine

end module