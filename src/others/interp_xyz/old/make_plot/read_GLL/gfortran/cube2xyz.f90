subroutine cube2xyz( &
  xyza, uvw, & 
  xyz, DuvwDxyz )
!input
! xyza: anchor points
! uvw: local coordinate 
!output
! xyz: map uvw to physical space
! DuvwDxyz

  use constants,only: NDIM,NGNOD
  use user_par,only: RK

  implicit none

  ! subroutine arguments
  real(RK), dimension(NDIM,NGNOD) :: xyza
  real(RK), dimension(NDIM) :: uvw

  real(RK), dimension(NDIM) :: xyz 
  real(RK), dimension(NDIM,NDIM) :: DuvwDxyz

  ! local variables
  real(RK), dimension(NDIM) :: lag1, lag2, lag3 
  real(RK), dimension(NDIM) :: lag1p, lag2p, lag3p
  real(RK), dimension(NGNOD) :: shape3D 
  real(RK), dimension(NGNOD,NDIM) :: dershape3D
  real(RK) :: jacobian
  real(RK), dimension(NDIM,NDIM) :: DxyzDuvw

  ! lagrange polynomials of order 3 on [-1,1], with collocation points: -1,0,1 
  lag1 = uvw*(uvw-1.0)/2.0
  lag2 = 1.0-uvw**2
  lag3 = uvw*(uvw+1.0)/2.0
  
  ! derivative of lagrange polynomials
  lag1p = uvw-0.5
  lag2p = -2.0*uvw
  lag3p = uvw+0.5
  
  ! construct the shape function
  shape3D = (/ &
       ! corner center
       lag1(1)*lag1(2)*lag1(3), &
       lag3(1)*lag1(2)*lag1(3), &
       lag3(1)*lag3(2)*lag1(3), & 
       lag1(1)*lag3(2)*lag1(3), &
       lag1(1)*lag1(2)*lag3(3), &
       lag3(1)*lag1(2)*lag3(3), &
       lag3(1)*lag3(2)*lag3(3), &
       lag1(1)*lag3(2)*lag3(3), &
       ! edge center
       lag2(1)*lag1(2)*lag1(3), &
       lag3(1)*lag2(2)*lag1(3), &
       lag2(1)*lag3(2)*lag1(3), &
       lag1(1)*lag2(2)*lag1(3), &
       lag1(1)*lag1(2)*lag2(3), &
       lag3(1)*lag1(2)*lag2(3), &
       lag3(1)*lag3(2)*lag2(3), &
       lag1(1)*lag3(2)*lag2(3), & 
       lag2(1)*lag1(2)*lag3(3), &
       lag3(1)*lag2(2)*lag3(3), &
       lag2(1)*lag3(2)*lag3(3), &
       lag1(1)*lag2(2)*lag3(3), &
       ! face center
       lag2(1)*lag2(2)*lag1(3), &
       lag2(1)*lag1(2)*lag2(3), &
       lag3(1)*lag2(2)*lag2(3), &
       lag2(1)*lag3(2)*lag2(3), &
       lag1(1)*lag2(2)*lag2(3), &
       lag2(1)*lag2(2)*lag3(3), &
       ! body center
       lag2(1)*lag2(2)*lag2(3) /)
                            
  ! derivative of the shape function
  ! corner center
  dershape3D( 1,:) = (/ lag1p(1)*lag1(2)*lag1(3), lag1(1)*lag1p(2)*lag1(3), lag1(1)*lag1(2)*lag1p(3) /)
  dershape3D( 2,:) = (/ lag3p(1)*lag1(2)*lag1(3), lag3(1)*lag1p(2)*lag1(3), lag3(1)*lag1(2)*lag1p(3) /)
  dershape3D( 3,:) = (/ lag3p(1)*lag3(2)*lag1(3), lag3(1)*lag3p(2)*lag1(3), lag3(1)*lag3(2)*lag1p(3) /) 
  dershape3D( 4,:) = (/ lag1p(1)*lag3(2)*lag1(3), lag1(1)*lag3p(2)*lag1(3), lag1(1)*lag3(2)*lag1p(3) /)
  dershape3D( 5,:) = (/ lag1p(1)*lag1(2)*lag3(3), lag1(1)*lag1p(2)*lag3(3), lag1(1)*lag1(2)*lag3p(3) /)
  dershape3D( 6,:) = (/ lag3p(1)*lag1(2)*lag3(3), lag3(1)*lag1p(2)*lag3(3), lag3(1)*lag1(2)*lag3p(3) /)
  dershape3D( 7,:) = (/ lag3p(1)*lag3(2)*lag3(3), lag3(1)*lag3p(2)*lag3(3), lag3(1)*lag3(2)*lag3p(3) /)
  dershape3D( 8,:) = (/ lag1p(1)*lag3(2)*lag3(3), lag1(1)*lag3p(2)*lag3(3), lag1(1)*lag3(2)*lag3p(3) /)
  ! edge center
  dershape3D( 9,:) = (/ lag2p(1)*lag1(2)*lag1(3), lag2(1)*lag1p(2)*lag1(3), lag2(1)*lag1(2)*lag1p(3) /)
  dershape3D(10,:) = (/ lag3p(1)*lag2(2)*lag1(3), lag3(1)*lag2p(2)*lag1(3), lag3(1)*lag2(2)*lag1p(3) /)
  dershape3D(11,:) = (/ lag2p(1)*lag3(2)*lag1(3), lag2(1)*lag3p(2)*lag1(3), lag2(1)*lag3(2)*lag1p(3) /)
  dershape3D(12,:) = (/ lag1p(1)*lag2(2)*lag1(3), lag1(1)*lag2p(2)*lag1(3), lag1(1)*lag2(2)*lag1p(3) /)
  dershape3D(13,:) = (/ lag1p(1)*lag1(2)*lag2(3), lag1(1)*lag1p(2)*lag2(3), lag1(1)*lag1(2)*lag2p(3) /)
  dershape3D(14,:) = (/ lag3p(1)*lag1(2)*lag2(3), lag3(1)*lag1p(2)*lag2(3), lag3(1)*lag1(2)*lag2p(3) /)
  dershape3D(15,:) = (/ lag3p(1)*lag3(2)*lag2(3), lag3(1)*lag3p(2)*lag2(3), lag3(1)*lag3(2)*lag2p(3) /)
  dershape3D(16,:) = (/ lag1p(1)*lag3(2)*lag2(3), lag1(1)*lag3p(2)*lag2(3), lag1(1)*lag3(2)*lag2p(3) /) 
  dershape3D(17,:) = (/ lag2p(1)*lag1(2)*lag3(3), lag2(1)*lag1p(2)*lag3(3), lag2(1)*lag1(2)*lag3p(3) /)
  dershape3D(18,:) = (/ lag3p(1)*lag2(2)*lag3(3), lag3(1)*lag2p(2)*lag3(3), lag3(1)*lag2(2)*lag3p(3) /)
  dershape3D(19,:) = (/ lag2p(1)*lag3(2)*lag3(3), lag2(1)*lag3p(2)*lag3(3), lag2(1)*lag3(2)*lag3p(3) /)
  dershape3D(20,:) = (/ lag1p(1)*lag2(2)*lag3(3), lag1(1)*lag2p(2)*lag3(3), lag1(1)*lag2(2)*lag3p(3) /)
  ! face center
  dershape3D(21,:) = (/ lag2p(1)*lag2(2)*lag1(3), lag2(1)*lag2p(2)*lag1(3), lag2(1)*lag2(2)*lag1p(3) /)
  dershape3D(22,:) = (/ lag2p(1)*lag1(2)*lag2(3), lag2(1)*lag1p(2)*lag2(3), lag2(1)*lag1(2)*lag2p(3) /)
  dershape3D(23,:) = (/ lag3p(1)*lag2(2)*lag2(3), lag3(1)*lag2p(2)*lag2(3), lag3(1)*lag2(2)*lag2p(3) /)
  dershape3D(24,:) = (/ lag2p(1)*lag3(2)*lag2(3), lag2(1)*lag3p(2)*lag2(3), lag2(1)*lag3(2)*lag2p(3) /)
  dershape3D(25,:) = (/ lag1p(1)*lag2(2)*lag2(3), lag1(1)*lag2p(2)*lag2(3), lag1(1)*lag2(2)*lag2p(3) /)
  dershape3D(26,:) = (/ lag2p(1)*lag2(2)*lag3(3), lag2(1)*lag2p(2)*lag3(3), lag2(1)*lag2(2)*lag3p(3) /)
  ! body center
  dershape3D(27,:) = (/ lag2p(1)*lag2(2)*lag2(3), lag2(1)*lag2p(2)*lag2(3), lag2(1)*lag2(2)*lag2p(3) /)

  ! xyz and Dxyz/Duvw
  xyz = matmul(xyza,shape3D)
  DxyzDuvw = matmul(xyza,dershape3D)

  ! adj(Dxyz/Duvw) -> DuvwDxyz
  DuvwDxyz(1,1) =   DxyzDuvw(2,2)*DxyzDuvw(3,3)-DxyzDuvw(3,2)*DxyzDuvw(2,3)
  DuvwDxyz(2,1) = -(DxyzDuvw(2,1)*DxyzDuvw(3,3)-DxyzDuvw(3,1)*DxyzDuvw(2,3))
  DuvwDxyz(3,1) =   DxyzDuvw(2,1)*DxyzDuvw(3,2)-DxyzDuvw(3,1)*DxyzDuvw(2,2)
 
  DuvwDxyz(1,2) = -(DxyzDuvw(1,2)*DxyzDuvw(3,3)-DxyzDuvw(3,2)*DxyzDuvw(1,3))
  DuvwDxyz(2,2) =   DxyzDuvw(1,1)*DxyzDuvw(3,3)-DxyzDuvw(3,1)*DxyzDuvw(1,3)
  DuvwDxyz(3,2) = -(DxyzDuvw(1,1)*DxyzDuvw(3,2)-DxyzDuvw(3,1)*DxyzDuvw(1,2))
  
  DuvwDxyz(1,3) =   DxyzDuvw(1,2)*DxyzDuvw(2,3)-DxyzDuvw(2,2)*DxyzDuvw(1,3)
  DuvwDxyz(2,3) = -(DxyzDuvw(1,1)*DxyzDuvw(2,3)-DxyzDuvw(2,1)*DxyzDuvw(1,3))
  DuvwDxyz(3,3) =   DxyzDuvw(1,1)*DxyzDuvw(2,2)-DxyzDuvw(2,1)*DxyzDuvw(1,2)

  ! jacobian = det(Dxyz/Duvw)
  jacobian = DxyzDuvw(1,1)*DuvwDxyz(1,1) &
           + DxyzDuvw(1,2)*DuvwDxyz(2,1) &
           + DxyzDuvw(1,3)*DuvwDxyz(3,1)
  if (jacobian<=0.0) stop '3D Jacobian undefined'

  ! Duvw/Dxyz = inv(Dxyz/Duvw) = adj(DxyzDuvw)/det(DxyzDuvw)  
  DuvwDxyz = DuvwDxyz/jacobian

end subroutine cube2xyz