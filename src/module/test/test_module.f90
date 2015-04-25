module constants

implicit none

integer :: nx

contains

subroutine set_dimension(n)

  integer, intent(in) :: n

  nx = n

end subroutine

end module constants 

!-------------------

module application

use constants

implicit none

contains

subroutine set_value(x)

  real, intent(inout) :: x(nx)
  
  ! local variables
  ! interesting
  integer :: ny
  real :: y(nx)
  integer :: z(5)

  ny = size(y)

  print *, "moduel:application: ny=", ny

  y = 1.0
  z = -12345

  x = y

end subroutine

end module application

!-------------------

program test

  use constants
  use application

  implicit none

  real, allocatable :: x(:)
  integer :: n_in

  character(len=10) :: args

  call get_command_argument(1, args)
  if (trim(args) == '') then
    args = "5"
  endif
  read(args,*) n_in

  call set_dimension(n_in)

  print *, "nx=", nx

  allocate(x(nx))

  call set_value(x)

  print *, "x=", x
  print *, "x(nx+1)=", x(nx+1)

end program test

! compile this as gfortran -g -fbounds-check ?.f90
