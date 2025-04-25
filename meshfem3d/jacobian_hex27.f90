!///////////////////////////////////////////////////////////////////////////////

subroutine xyz2cube_bounded_hex27(xyz_anchor, xyz, uvw, misloc, flag_inside)
!-mapping a given point in physical space (xyz) to the 
! reference cube (uvw) for a 27-node hexahedron element (8 vertex, 12 egde
! centers, 6 face centers and 1 body center),
! and also flag whether the point is inside the cube
! if the point lies outside the element, calculate the bounded (xi,eta,gamma)
! inside or on the surface of the reference unit cube.
!
!-inputs:
! (real) xyz_anchor(3,27): anchor points of the element
! (real) xyz(3): coordinates of the target point
!
!-outputs:
! (real) uvw(3): local coordinates in reference cube
! (real) misloc: location misfit abs(xyz - XYZ(uvw))
! (logical) flag_inside: flag whether the target point locates inside the element

  implicit none

  ! input/output
  real(kind=8), intent(in) :: xyz_anchor(3,27)
  real(kind=8), intent(in) :: xyz(3)

  real(kind=8), intent(out) :: uvw(3)
  real(kind=8), intent(out) :: misloc
  logical, intent(out) :: flag_inside

  ! local variables

  ! number of iterations used to locate point inside one element 
  integer, parameter :: niter = 5

  integer :: iter
  real(kind=8), dimension(3) :: xyzi ! iteratively improved xyz
  real(kind=8), dimension(3,3) :: DuvwDxyz
  real(kind=8), dimension(3) :: dxyz, duvw
  real(kind=8), parameter ::  ZERO=0.d0, ONE=1.d0, MINUS_ONE=-1.d0

  ! initialize 
  uvw = ZERO
  flag_inside = .true.

  ! iteratively update local coordinate uvw to approach the target xyz
  do iter = 1, niter

    ! predicted xyzi and Jacobian for the current uvw
    call jacobian_hex27(xyz_anchor, uvw, xyzi, DuvwDxyz)

    ! compute difference
    dxyz = xyz - xyzi

    ! compute increments
    duvw = matmul(DuvwDxyz, dxyz)

    ! update values
    uvw = uvw + duvw

    ! limit inside the cube
    if (any(uvw < MINUS_ONE .or. uvw > ONE)) then 
      where (uvw < MINUS_ONE) uvw = MINUS_ONE
      where (uvw > ONE) uvw = ONE
      ! set is_inside to false based on the last iteration
      if (iter == niter) then
        flag_inside = .false.
      endif
    endif

  enddo ! do iter_loop = 1,NUM_ITER
  
  ! calculate the predicted position 
  call jacobian_hex27(xyz_anchor, uvw, xyzi, DuvwDxyz)

  ! residual distance from the target point
  misloc = sqrt(sum((xyz-xyzi)**2))

end subroutine xyz2cube_bounded_hex27

!///////////////////////////////////////////////////////////////////////////////

subroutine jacobian_hex27(xyz_anchor, uvw, xyz, DuvwDxyz)
! compute 3D jacobian at a given point for a 27-node element
! map from local coordinate (uvw) to physical position (xyz)
! the shape the element is defined by the anchor points (xyz_anchor)
!
!-input
! xyz_anchor(3,27): xyz of anchor points, the order of the 27 nodes must be from
!   subroutine anchor_index_hex27()
! uvw(3): local coordinate 
!
!-output
! xyz(3): map uvw to physical space
! DuvwDxyz(3,3): jacobian matrix

  implicit none

  ! input/output
  real(kind=8), intent(in) :: xyz_anchor(3,27), uvw(3)
  real(kind=8), intent(out) :: xyz(3), DuvwDxyz(3,3)
 
  ! local variables
  real(kind=8), dimension(3) :: lag1, lag2, lag3 
  real(kind=8), dimension(3) :: lag1p, lag2p, lag3p
  real(kind=8), dimension(27) :: shape3D 
  real(kind=8), dimension(27,3) :: dershape3D
  real(kind=8) :: jacobian
  real(kind=8), dimension(3,3) :: DxyzDuvw

  ! lagrange polynomials of order 3 on [-1,1], with collocation points: -1,0,1 
  lag1 = uvw * (uvw - 1.d0) / 2.d0
  lag2 = 1.d0 - uvw**2
  lag3 = uvw * (uvw + 1.d0) / 2.d0
  
  ! derivative of lagrange polynomials
  lag1p = uvw - 0.5d0
  lag2p = -2.d0 * uvw
  lag3p = uvw + 0.5d0
  
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
  xyz = matmul(xyz_anchor, shape3D)
  DxyzDuvw = matmul(xyz_anchor, dershape3D)

  ! adjoint matrix: adj(Dxyz/Duvw)
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
  jacobian =  DxyzDuvw(1,1)*DuvwDxyz(1,1) &
            + DxyzDuvw(1,2)*DuvwDxyz(2,1) &
            + DxyzDuvw(1,3)*DuvwDxyz(3,1)

  if (jacobian <= 0.d0) then
    print *, "[ERROR] jacobian_hex27: 3D Jacobian undefined jacobian=", jacobian
    print *, "xyz_anchor(3,27)=", xyz_anchor
    print *, "uvw(3)=", uvw
    stop
  endif

  ! inverse matrix: Duvw/Dxyz = inv(Dxyz/Duvw) = adj(DxyzDuvw)/det(DxyzDuvw)  
  DuvwDxyz = DuvwDxyz / jacobian
 
end subroutine jacobian_hex27


subroutine anchor_index_hex27(NGLLX,NGLLY,NGLLZ,iax,iay,iaz)
!-return indices of hex27 anchor points in a GLL element used for jacobian_hex27
!
!-input: 
!   NGLLX,Y,Z: number of GLL points along x,y,z axes
!
!-output:
!   iax,iay,iaz(27): indices on x,y,z axes of the 27 nodes

  implicit none

  ! input/output
  integer, intent(in) :: NGLLX, NGLLY, NGLLZ
  integer, intent(out) :: iax(27), iay(27), iaz(27)

  ! local vars
  integer :: MIDX, MIDY, MIDZ
  integer :: ENDX, ENDY, ENDZ

  ! mid-point index
  MIDX = (NGLLX-1)/2
  MIDY = (NGLLY-1)/2
  MIDZ = (NGLLZ-1)/2
  ! end-point index
  ENDX = NGLLX - 1
  ENDY = NGLLY - 1
  ENDZ = NGLLZ - 1

  ! corner nodes
  iax(1) = 0;    iay(1) = 0;    iaz(1) = 0
  iax(2) = ENDX; iay(2) = 0;    iaz(2) = 0
  iax(3) = ENDX; iay(3) = ENDY; iaz(3) = 0
  iax(4) = 0;    iay(4) = ENDY; iaz(4) = 0
  iax(5) = 0;    iay(5) = 0;    iaz(5) = ENDZ
  iax(6) = ENDX; iay(6) = 0;    iaz(6) = ENDZ
  iax(7) = ENDX; iay(7) = ENDY; iaz(7) = ENDZ
  iax(8) = 0;    iay(8) = ENDY; iaz(8) = ENDZ

  ! midside nodes (nodes located in the middle of an edge)
  iax( 9) = MIDX; iay( 9) = 0;    iaz( 9) = 0
  iax(10) = ENDX; iay(10) = MIDY; iaz(10) = 0
  iax(11) = MIDX; iay(11) = ENDY; iaz(11) = 0
  iax(12) = 0;    iay(12) = MIDY; iaz(12) = 0
  iax(13) = 0;    iay(13) = 0;    iaz(13) = MIDZ
  iax(14) = ENDX; iay(14) = 0;    iaz(14) = MIDZ
  iax(15) = ENDX; iay(15) = ENDY; iaz(15) = MIDZ
  iax(16) = 0;    iay(16) = ENDY; iaz(16) = MIDZ
  iax(17) = MIDX; iay(17) = 0;    iaz(17) = ENDZ
  iax(18) = ENDX; iay(18) = MIDY; iaz(18) = ENDZ
  iax(19) = MIDX; iay(19) = ENDY; iaz(19) = ENDZ
  iax(20) = 0;    iay(20) = MIDY; iaz(20) = ENDZ

  ! side center nodes (nodes located in the middle of a face)
  iax(21) = MIDX; iay(21) = MIDY; iaz(21) = 0
  iax(22) = MIDX; iay(22) = 0;    iaz(22) = MIDZ
  iax(23) = ENDX; iay(23) = MIDY; iaz(23) = MIDZ
  iax(24) = MIDX; iay(24) = ENDY; iaz(24) = MIDZ
  iax(25) = 0;    iay(25) = MIDY; iaz(25) = MIDZ
  iax(26) = MIDX; iay(26) = MIDY; iaz(26) = ENDZ

  ! center node (barycenter of the eight corners)
  iax(27) = MIDX; iay(27) = MIDY; iaz(27) = MIDZ

end subroutine anchor_index_hex27
