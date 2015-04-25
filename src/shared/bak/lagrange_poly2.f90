!interpolation coefficients: L_i(x)
subroutine lagrange_poly(x, L)

  use constants, only: XGLL, NGLLX, RK=>SIZE_REAL

  implicit none

  ! subroutine arguments
  real(RK) :: x
  real(RK), dimension(NGLLX) :: L

  ! local variables
  integer :: i
  integer, dimension(NGLLX) :: ind
  real(RK), dimension(NGLLX) :: xx, yy

  ! L_i(x)
  ind = (/(i,i=1,NGLLX)/)
  xx = x-XGLL
  do i = 1,NGLLX
    yy = XGLL(i)-XGLL
    L(i) = product(xx/yy,MASK=(ind/=i))
  end do

end subroutine lagrange_poly