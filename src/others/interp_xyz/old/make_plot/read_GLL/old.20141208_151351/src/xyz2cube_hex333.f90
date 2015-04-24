! invert map a given point in physical space (xyz) to the 
! 3x3x3 reference cube (xi,eta,gamma), 
! and also flag if the point is inside the cube

subroutine xyz2cube_hex333 &
( xa,ya,za, x,y,z, &
  xi,eta,gamma, dist_res, is_inside )

	implicit none

	include "constants.h"	
	include "values_from_mesher.h"	
	
	!==========================================================
	! specification of subroutine arguments 
	!==========================================================
	
	! target point coordinate
	double precision, intent(in) :: x, y, z
	
	! anchor points of the 3x3x3 hex 
	real(kind=CUSTOM_REAL), dimension(NGNOD) :: xa, ya, za
	
	! initial coordinate in the reference cube
	double precision, intent(inout) :: xi, eta, gamma
	
	! flag if point is inside the cube
	double precision, intent(out) :: dist_res
	integer, intent(out) :: is_inside
	
	!==========================================================
	! define local variables
	!==========================================================
	
	integer :: iter_loop
	double precision :: xcur,ycur,zcur ! iteratively improve current xyz

	double precision :: jacobian,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
	double precision :: dx,dy,dz,dxi,deta,dgamma
	
	!==========================================================
	! find the cube coordinate 
	!==========================================================

	is_inside = 1
	
	! iteratively update (xi,eta,gamma) to approach the current x,y,z to the target xyz
	do iter_loop = 1,NUM_ITER
	
	  ! compute xyz and Jacobian for the current point (xi,eta,gamma)
	  call cube2xyz_hex333(xa,ya,za,xi,eta,gamma, &
	                       xcur,ycur,zcur, &
	                       jacobian,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

		! if jacobian <= 0 exit
		if (jacobian <= ZERO ) then
			is_inside = 0	
			return
		endif 
		
	  ! compute difference
	  dx = (x - xcur)
	  dy = (y - ycur)
	  dz = (z - zcur)
	
	  ! compute increments
	  dxi  = xix*dx + xiy*dy + xiz*dz
	  deta = etax*dx + etay*dy + etaz*dz
	  dgamma = gammax*dx + gammay*dy + gammaz*dz
	
	  ! update values
	  xi = xi + dxi
	  eta = eta + deta
	  gamma = gamma + dgamma
	
	enddo ! do iter_loop = 1,NUM_ITER
	
	! final update the position
	call cube2xyz_hex333(xa,ya,za,xi,eta,gamma, &
	                      xcur,ycur,zcur, &
	                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
	
	! residual distance from the target point
	dist_res = dsqrt((x-xcur)**2+(y-ycur)**2+(z-zcur)**2)

	! flag if point is inside the cube
	if (xi<-1.d0 .or. xi>1.d0 .or. & 
	    eta<-1.d0 .or. eta>1.d0 .or. &
	    gamma<-1.d0 .or. gamma>1.d0 ) then
	  is_inside = 0
	endif

end subroutine xyz2cube_hex333
