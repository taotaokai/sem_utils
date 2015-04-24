! invert mapping a given point in physical space (xyz) to the 
! 3x3x3 reference cube (xi,eta,gamma), 
! and also flag if the point is inside the cube
! if the point lies outside the element, calculate the bounded (xi,eta,gamma)
! inside or on the surface of the reference unit cube.

subroutine xyz2cube_bounded ( &
  xyza, xyz, &
  uvw, res, inside )
!input
! xyza: shape function anchor points
! xyz: coordinates of the target point
!output
! uvw: local coordinates in reference cube, defined by anchor points xyza
! res: residual distance of target and predicted points from (uvw)
! inside: flag whether the target point is inside the host element

  use constants,only: NUM_ITER, NGNOD, NDIM
  use user_par,only: RK

  implicit none

  !==========================================================
  ! specification of subroutine arguments 
  !==========================================================
  ! anchor points of the 3x3x3 hex
  real(RK), dimension(NDIM,NGNOD), intent(in) :: xyza

  ! target point coordinate
  real(RK), dimension(NDIM), intent(in) :: xyz
  
  ! local coordinate in the reference cube
  real(RK), dimension(NDIM) :: uvw
  
  ! flag if point is inside the cube
  real(RK) :: res
  logical :: inside
  
  !==========================================================
  ! define local variables
  !==========================================================
  
  integer :: iter_loop
  real(RK), dimension(NDIM) :: xyzi ! iteratively improved xyz
  real(RK), dimension(NDIM,NDIM) :: DuvwDxyz
  real(RK), dimension(NDIM) :: dxyz,duvw

  !==========================================================
  ! find the cube coordinate 
  !==========================================================
  uvw = 0.0
  inside = .true.
  ! iteratively update local coordinate uvw to approach the target xyz
  do iter_loop = 1,NUM_ITER
    ! predicted xyzi and Jacobian for the current uvw
    call cube2xyz(xyza,uvw, xyzi,DuvwDxyz)
    ! compute difference
    dxyz = xyz-xyzi
    ! compute increments
    duvw = matmul(DuvwDxyz,dxyz)
    ! update values
    uvw = uvw + duvw
    ! limit inside the cube
    if (any(uvw<=-1.0 .or. uvw>=1.0)) then 
      where (uvw<-1.0) uvw = -1.0
      where (uvw>1.0) uvw = 1.0
      ! flag inside in the last iteration step
      if (iter_loop==NUM_ITER) inside=.false.
    end if
  end do ! do iter_loop = 1,NUM_ITER
  
  ! calculate the predicted position 
  call cube2xyz(xyza,uvw,xyzi,DuvwDxyz)

  ! residual distance from the target point
  res = sqrt(sum((xyz-xyzi)**2))

end subroutine xyz2cube_bounded
