! subroutine to compute the Lagrange interpolants based upon the GLL points
! at any point xi in [-1,1]

subroutine lagrange(xi,NGLL,xigll,lxi)

  implicit none
  
  double precision, intent(in) :: xi
  integer, intent(in) :: NGLL
  double precision, dimension(NGLL), intent(in) :: xigll
  double precision, dimension(NGLL), intent(out) :: lxi
  
  integer :: dgr,i
  double precision :: prod1,prod2
  
  do dgr = 1,NGLL
     prod1 = 1.0d0
     prod2 = 1.0d0
  
     do i = 1,NGLL
        if (i /= dgr) then
           prod1 = prod1 * (xi         - xigll(i))
           prod2 = prod2 * (xigll(dgr) - xigll(i))
        endif
     enddo
  
     lxi(dgr) = prod1 / prod2
  enddo
  
end subroutine lagrange
