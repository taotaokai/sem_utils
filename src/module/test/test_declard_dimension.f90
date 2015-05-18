! test declaring array dimensions using an integer array of dims

program test

  implicit none

  integer :: i, j

  !integer, dimension(2), parameter :: dims = 2

  integer, dimension(2, 2) :: array

  do i = 1, 2
    do j = 1, 2
      array(i,j) = i+j
    enddo
  enddo

  print *, "array(:,:)=", array

end 