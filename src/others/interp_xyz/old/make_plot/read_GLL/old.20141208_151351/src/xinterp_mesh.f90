! locate point in the cube

program xinterp_mesh

	implicit none

	include "constants.h"
	include "values_from_mesher.h"

	!use constants,only: CUSTOM_REAL, NGLOB_CRUST_MANTLE, GAUSSALPHA, GAUSSBETA

  !implicit none

  !==========================================================
  ! user parameters 
  !==========================================================

  ! element control points
	!integer, parameter :: NGLLX=5!, NGNOD=27 ! must be this value

	! interpolation points
	integer, parameter :: MAXNUMPOINT=20000

	!==========================================================
  ! local variables 
  !==========================================================

  ! database directory and model name
  character(len=256) :: databasedir, modelname

  ! region id and process number
  integer :: iregion, nproc

	! others
	integer :: iproc, ispec, ipoint, iglob
	integer :: ia, i,j,k, ierr
	double precision :: xc,yc,zc, x,y,z, dist, xi,eta,gamma, dist_residual

	! mesh data
  integer :: nspec, nglob
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xmesh, ymesh, zmesh
  integer, dimension(NGLLX,NGLLX,NGLLX,NSPEC_CRUST_MANTLE) :: ibool

  ! model data
  real(kind=CUSTOM_REAL) :: vmesh(NGLLX,NGLLX,NGLLX,NSPEC_CRUST_MANTLE)

  ! initial cube coordinate
  double precision :: xi0, eta0, gamma0
  integer :: is_inside

	! gll
  double precision, dimension(NGLLX) :: xigll, wxgll, lxi, leta, lgamma

	! distance threshold
	double precision :: dist_threshold

  ! element control points
  real(kind=CUSTOM_REAL) :: xa(NGNOD), ya(NGNOD), za(NGNOD) 
  integer, dimension(NGNOD) :: ixa, iya, iza ! index of anchor points on hex333 in spec555

	! interpolation points
	integer :: npoint	
  real(kind=CUSTOM_REAL), dimension(MAXNUMPOINT) :: xpoint, ypoint, zpoint
	double precision :: vpoint
	character(len=30) :: FMTOUT

  !==========================================================
  ! do the real work 
  !==========================================================

	! read in program arguments
	call get_program_args(databasedir,modelname,nproc,iregion)

  ! read in point coordinates
	i = 1
	do
		read(*,*,iostat=ierr) xpoint(i),ypoint(i),zpoint(i)
		if (ierr /= 0) exit
		!write(*,*) "point read in: ",xpoint(i),ypoint(i),zpoint(i)
		i = i+1
	enddo
	npoint = i-1
 	!write(*,*) "npoint=",npoint

	! for 5x5x5 specfem, get the anchor/contrl point index (3x3x3 hex)
	call index_hex333_in_spec555(ixa,iya,iza)

	! (xi,eta,gamma) at the center of spec555
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
	xi0 = xigll((NGLLX+1)/2)
	eta0 = xigll((NGLLX+1)/2)
	gamma0 = xigll((NGLLX+1)/2);

	! output format
	FMTOUT = "(I5, 3F12.8, F12.8)"
  
	! loop each mesh part
	do iproc = 1,nproc
		
  	! read mesh data 
  	call read_mesh(iproc-1, iregion, databasedir, modelname, &
   	               nspec, nglob, xmesh, ymesh, zmesh, vmesh, ibool)

		!write(*,*) "nspec=",nspec

	  ! loop each element
	  do ispec = 1,nspec

			!write(*,*) "ispec=",ispec

			! get anchor point locations
			!write(*,*) "anchor points"
			do ia = 1,NGNOD
				iglob = ibool(ixa(ia),iya(ia),iza(ia),ispec);
				xa(ia) = xmesh(iglob)
				ya(ia) = ymesh(iglob)
				za(ia) = zmesh(iglob)
				!write(*,*) xa(ia),ya(ia),za(ia),vmesh(ixa(ia),iya(ia),iza(ia),ispec)
			enddo

			! element center
			xc = xa(27); yc = ya(27); zc = za(27)
			!write(*,"(F12.8,F12.8,F12.8)") xc,yc,zc

      ! distance threshold
			dist_threshold = 0
			do ia = 1,NGNOD
  			dist = dsqrt((xa(ia)-xc)**2+(ya(ia)-yc)**2+(za(ia)-zc)**2)
				if (dist > dist_threshold) dist_threshold = dist
			enddo

	    ! loop each point
			do ipoint = 1,npoint

				!write(*,*) "ipoint=",ipoint

				x = xpoint(ipoint); y = ypoint(ipoint); z = zpoint(ipoint);

				! first exclude points too far away from element center
  			dist = dsqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
				!write(*,*) "----- dist from center: ",dist, dist_threshold
			  if (dist > dist_threshold) cycle

				!write(*,*) "===== passed first check"

			  ! second, compute the cube coordinate

			  xi = xi0; eta = eta0; gamma = gamma0 ! set initial point in the center
			  call xyz2cube_hex333(xa,ya,za,x,y,z,xi,eta,gamma,dist_residual,is_inside)
				!write(*,*) "----- xi/eta/gamma= ",xi,eta,gamma,dist_residual
				if (is_inside /= 1) cycle 

				! interpolate point value
				call lagrange(xi, NGLLX, xigll, lxi)
				call lagrange(eta, NGLLX, xigll, leta)
				call lagrange(gamma, NGLLX, xigll, lgamma)
				vpoint = 0.d0
				do k = 1,NGLLX
				  do j = 1,NGLLX
				    do i = 1,NGLLX
				      vpoint = vpoint + lxi(i)*leta(j)*lgamma(k)*vmesh(i,j,k,ispec)
				    enddo
				  enddo
				enddo

  			! write out model values for each point
				write(*,FMTOUT) ipoint,x,y,z,vpoint

			enddo ! do ipoint
		enddo ! do ispec
	enddo ! do iproc
	
end program xinterp_mesh
