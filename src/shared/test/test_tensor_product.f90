program test

  integer :: a(3,1)
  integer :: b(1,3)
  integer :: c(3,3)

  a(1,1) = 1
  a(2,1) = 2
  a(3,1) = 3

  b(1,1) = 2
  b(1,2) = 3
  b(1,3) = 4

  c = matmul(a,b)

  do i = 1,3
    do j = 1,3
      print '(a)', i, j, c(i,j)
    enddo
  enddo

end program