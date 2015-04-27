module sem_mesh

!----------------------------------
! specification part
!----------------------------------

  use sem_constants
  use sem_io

  implicit none

  private

  ! for calculation of volume correspoding to each gll point
  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll
  ! integration weights on GLL points in the unit cube [-1,1]^NDIM
  real(kind=CUSTOM_REAL) :: wgll_cube(NGLLX,NGLLY,NGLLZ)

  ! mesh data type 
  type, public :: mesh
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyz
    integer, dimension(:,:,:,:), allocatable :: ibool
    integer, dimension(:), allocatable :: idoubling
    logical, dimension(:), allocatable :: ispec_is_tiso 
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  end type mesh

  ! public operations
  public :: sem_mesh_init
  public :: sem_mesh_read
  public :: sem_mesh_get_gll_volume
  public :: sem_mesh_locate_xyz

!----------------------------------
! implementation part
!----------------------------------
contains

!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_init(mesh_data)

  type (mesh), intent(inout) :: mesh_data

  if (allocated(mesh_data%xyz)) then
    deallocate(mesh_data%xyz, &
               mesh_data%ibool, &
               mesh_data%idoubling, &
               mesh_data%ispec_is_tiso, &
               mesh_data%xix, &
               mesh_data%xiy, &
               mesh_data%xiz, &
               mesh_data%etax, &
               mesh_data%etay, &
               mesh_data%etaz, &
               mesh_data%gammax, &
               mesh_data%gammay, &
               mesh_data%gammaz)
  end if

  allocate(mesh_data%xyz(1:3,1:NGLOB), &
           mesh_data%ibool(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%idoubling(1:NSPEC), &
           mesh_data%ispec_is_tiso(1:NSPEC), &
           mesh_data%xix(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%xiy(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%xiz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%etax(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%etay(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%etaz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%gammax(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%gammay(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%gammaz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC))

  ! compute local coordinates of gll points, integration weight
  call compute_wgll_cube()

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_read(basedir, iproc, iregion, mesh_data)

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion
  type (mesh), intent(inout) :: mesh_data

  integer :: ival

  call sem_io_open_proc_file_for_read(IIN,basedir,iproc,iregion,'solver_data')

  ! nspec
  read(IIN) ival
  if (ival /= NSPEC) then
    write(*,*) "ERROR:sem_mesh: inconsistent NSPEC value!"
    write(*,*) "- solver_data.bin: nspec=", ival
    write(*,*) "- sem_constants: NSPEC=", NSPEC
    stop
  end if

  ! nglob
  read(IIN) ival
  if (ival /= NGLOB) then
    write(*,*) "ERROR:sem_mesh: inconsistent NGLOB value!"
    write(*,*) "- solver_data.bin: nglob=", ival
    write(*,*) "- sem_constants: NGLOB=", NGLOB 
    stop
  end if

  read(IIN) mesh_data%xyz(1,1:NGLOB)
  read(IIN) mesh_data%xyz(2,1:NGLOB)
  read(IIN) mesh_data%xyz(3,1:NGLOB)

  read(IIN) mesh_data%ibool(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%idoubling(1:NSPEC)
  read(IIN) mesh_data%ispec_is_tiso(1:NSPEC)

  read(IIN) mesh_data%xix(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%xiy(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%xiz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%etax(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%etay(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%etaz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%gammax(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%gammay(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)
  read(IIN) mesh_data%gammaz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC)

  close(IIN)

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_get_gll_volume(mesh_data, gll_volume)

! compute volume corresponding to each gll point in each spectral element
! gll_volume = jacobian * wgll_cube
! jacobian = det(d(x,y,z)/d(xi,eta,gamma)), volumetric deformation ratio from cube to actual element 

  type (mesh), intent(in) :: mesh_data
  real(kind=CUSTOM_REAL), intent(out) :: gll_volume(NGLLX,NGLLY,NGLLZ,NSPEC)

  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: jacobianl

  ! calculate gll quad weight in the reference cube
  !call compute_wgll_cube()

  ! calculate volume integration weight on gll
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! gets derivatives of ux, uy and uz with respect to x, y and z
                xixl = mesh_data%xix(i,j,k,ispec)
                xiyl = mesh_data%xiy(i,j,k,ispec)
                xizl = mesh_data%xiz(i,j,k,ispec)
              etaxl = mesh_data%etax(i,j,k,ispec)
              etayl = mesh_data%etay(i,j,k,ispec)
              etazl = mesh_data%etaz(i,j,k,ispec)
          gammaxl = mesh_data%gammax(i,j,k,ispec)
          gammayl = mesh_data%gammay(i,j,k,ispec)
          gammazl = mesh_data%gammaz(i,j,k,ispec)

          ! jacobian: det( d(x,y,z)/d(xi,eta,gamma))
          jacobianl = 1.0_CUSTOM_REAL / &
            ( xixl*(etayl*gammazl-etazl*gammayl)  &
             -xiyl*(etaxl*gammazl-etazl*gammaxl)  &
             +xizl*(etaxl*gammayl-etayl*gammaxl))

          ! gll volume: jacobian * wgll_cube 
          gll_volume(i,j,k,ispec) = jacobianl * wgll_cube(i,j,k)

        enddo
      enddo
    enddo
  enddo
  
end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine compute_wgll_cube()

! compute the integration weight for gll quadrature in a reference cube
! i.e. volume corresponding to gll point in the reference cube

  integer :: i, j, k

  call zwgljd(xigll, wxgll, NGLLX, GAUSSALPHA, GAUSSBETA)
  call zwgljd(yigll, wygll, NGLLY, GAUSSALPHA, GAUSSBETA)
  call zwgljd(zigll, wzgll, NGLLZ, GAUSSALPHA, GAUSSBETA)

  do k=1, NGLLZ
    do j=1, NGLLY
      do i=1, NGLLX
        wgll_cube(i,j,k) = real(wxgll(i)*wygll(j)*wzgll(k), kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_locate_xyz(mesh_data, npoint, xyz, uvw, hlagrange, misloc, elem_ind, stat_loc)
!-locate xyz(3,npoint) inside the mesh
!
!-input
! type(mesh) mesh_data
! xyz(3, npoint): coordinates to interpolate
!
!-output
! uvw(3, npoint): local coordinate inside the located element
! hlagrange(NDIM, NGLLX, NGLLY, NGLLZ, npoint): GLL interpolant value on uvw
! misloc(npoint): mis-location
! elem_ind(npoint): element index contains the located point
! stat_loc(npoint): location status (1: inside, 0: close to the element, -1: far away)

  type (mesh), intent(in) :: mesh_data
  integer, intent(in) :: npoint
  ! points to locate
  real(kind=CUSTOM_REAL), intent(in) :: xyz(NDIM, npoint)

  ! local coordinates
  real(kind=CUSTOM_REAL), intent(out) :: uvw(NDIM, npoint)
  real(kind=CUSTOM_REAL), intent(out) :: hlagrange(NGLLX, NGLLY, NGLLZ, npoint)
  ! mis-location: abs(xyz - XYZ(uvw))
  real(kind=CUSTOM_REAL), intent(out) :: misloc(npoint)
  ! index of the element inside which the point is located
  integer, intent(out) :: elem_ind(npoint)
  ! location status: 1 inside, 0 close, -1 away from a typical distance
  integer, intent(out) :: stat_loc(npoint)

  ! local varaibles
  integer :: i, ia, nsel, ispec, iglob, ipoint, isel
  integer :: RANGE_1_NSPEC(NSPEC), ind_sel(NSPEC), ind_sort(NSPEC)
  real(kind=CUSTOM_REAL) :: typical_size, HUGEVAL, dist_sq(NSPEC)
  real(kind=CUSTOM_REAL) :: max_x, min_x, max_y, min_y, max_z, min_z
  real(kind=CUSTOM_REAL) :: xyz1(3), uvw1(3), misloc1
  logical :: idx_sel(NSPEC), is_inside
  integer :: igllx, iglly, igllz
  double precision :: lagx(NGLLX), lagy(NGLLY), lagz(NGLLZ)

  integer, dimension(NGNOD) :: iaddx, iaddy, iaddz ! index of anchor points
  real(kind=CUSTOM_REAL) :: center_xyz(3, NSPEC)
  real(kind=CUSTOM_REAL) :: anchor_xyz(3, NGNOD, NSPEC)

  ! initialize some parameters

  ! search range: 5 * typical element width at surface
  typical_size = 5 * real((ANGULAR_WIDTH_XI_IN_DEGREES_VAL &
                * DEGREES_TO_RADIANS / NEX_XI_VAL) * R_UNIT_SPHERE &
                , kind=CUSTOM_REAL)
  HUGEVAL = huge(1.0_CUSTOM_REAL)
  RANGE_1_NSPEC = (/(i, i=1,NSPEC)/)

  max_x = maxval(center_xyz(1,:)) + typical_size
  min_x = minval(center_xyz(1,:)) - typical_size
  max_y = maxval(center_xyz(2,:)) + typical_size
  min_y = minval(center_xyz(2,:)) - typical_size
  max_z = maxval(center_xyz(3,:)) + typical_size
  min_z = minval(center_xyz(3,:)) - typical_size

  ! get index of anchor points in the GLL element
  call hex_nodes(iaddx,iaddy,iaddz)

  ! get anchor and center points of each GLL element 
  do ispec = 1, NSPEC
    do ia = 1, NGNOD
      iglob = mesh_data%ibool(iaddx(ia),iaddy(ia),iaddz(ia),ispec)
      anchor_xyz(:,ia,ispec) = mesh_data%xyz(:,iglob)
    enddo
    ! the last anchor point is the element center
    center_xyz(:,ispec) = mesh_data%xyz(:,iglob)
  enddo

  ! locate each point
  do ipoint = 1, npoint

    ! target point
    xyz1 = xyz(:,ipoint)

    ! skip points based on coordinate ranges
    if (     xyz1(1)<min_x .or. xyz1(1)>max_x &
        .or. xyz1(2)<min_y .or. xyz1(2)>max_y &
        .or. xyz1(3)<min_z .or. xyz1(3)>max_z) then
      stat_loc(ipoint) = -1
      elem_ind(ipoint) = -1
      cycle
    end if

    ! select elements close to the target point based on box coordinate range
    idx_sel =      (center_xyz(1,:)>xyz1(1)-typical_size) &
             .and. (center_xyz(1,:)<xyz1(1)+typical_size) &
             .and. (center_xyz(2,:)>xyz1(2)-typical_size) &
             .and. (center_xyz(2,:)<xyz1(2)+typical_size) &
             .and. (center_xyz(3,:)>xyz1(3)-typical_size) &
             .and. (center_xyz(3,:)<xyz1(3)+typical_size)

    nsel = count(idx_sel)

    if (nsel == 0) then 
      stat_loc(ipoint) = -1
      elem_ind(ipoint) = -1
      cycle
    endif

    ind_sel(1:nsel) = pack(range_1_NSPEC, idx_sel)

    ! sort distance
    dist_sq(1:nsel) =  (center_xyz(1,ind_sel(1:nsel)) - xyz1(1))**2 & 
                     + (center_xyz(2,ind_sel(1:nsel)) - xyz1(2))**2 &
                     + (center_xyz(3,ind_sel(1:nsel)) - xyz1(3))**2

    call heap_sort(nsel, dist_sq, ind_sort)

    ! loop all selected elements to find the one containing the target xyz 
    misloc(ipoint) = HUGEVAL

    do isel = 1, nsel

      ispec = ind_sel(ind_sort(isel))

      call xyz2cube_bounded(anchor_xyz(:,:,ispec),xyz1,uvw1,misloc1,is_inside)

      if (is_inside) then ! record this point and done
        stat_loc(ipoint) = 1
        elem_ind(ipoint) = ispec
        misloc(ipoint) = misloc1
        uvw(:,ipoint) = uvw1
        exit
      else ! record the nearest location
        if (misloc1 < misloc(ipoint)) then
          stat_loc(ipoint) = 0
          elem_ind(ipoint) = ispec
          misloc(ipoint) = misloc1
          uvw(:,ipoint) = uvw1
        endif
      endif !inside

    enddo ! isel

    ! set GLL interpolating values if located
    if (stat_loc(ipoint) /= -1) then
      call lagrange_poly(DBLE(uvw(1,ipoint)), NGLLX, xigll, lagx)
      call lagrange_poly(DBLE(uvw(2,ipoint)), NGLLY, yigll, lagy)
      call lagrange_poly(DBLE(uvw(3,ipoint)), NGLLZ, zigll, lagz)

      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            hlagrange(igllx,iglly,igllz,ipoint) &
              = real(lagx(igllx) * lagy(iglly) * lagz(igllz), kind=CUSTOM_REAL)
          enddo
        enddo
      enddo 

    endif

  enddo ! ipoint

end subroutine

!///////////////////////////////////////////////////////////////////////////////
subroutine lagrange_poly(x, ngll, xgll, lagrange)
!-get lagrange interpolation coefficients: L_i(x)
!
!-input
! x: coordinates of interpolating point
! ngll: number of colocation points
! xgll(ngll): coordinates of colocation points
!
!-output
! lagrange(ngll): interpolation coeff.
!

  double precision, intent(in) :: x
  integer, intent(in) :: ngll
  double precision, intent(in) :: xgll(ngll)

  double precision, intent(out) :: lagrange(ngll)

  ! local variables
  integer :: i
  integer, dimension(ngll) :: ind
  double precision, dimension(ngll) :: xx, yy

  ! lagrange(ngll) 
  ind = (/(i, i=1,ngll)/)
  xx = x - xgll

  do i = 1, ngll
    yy = xgll(i) - xgll
    lagrange(i) = product(xx/yy, mask=(ind/=i))
  end do

end subroutine lagrange_poly

!///////////////////////////////////////////////////////////////////////////////
subroutine xyz2cube_bounded(anchor_xyz, xyz, uvw, misloc, is_inside )
!-invert mapping a given point in physical space (xyz) to the 
! 3x3x3 reference cube (xi,eta,gamma), 
! and also flag if the point is inside the cube
! if the point lies outside the element, calculate the bounded (xi,eta,gamma)
! inside or on the surface of the reference unit cube.
!
!-input
! anchor_xyz: anchor points of the element
! xyz: coordinates of the target point
!
!-output
! uvw: local coordinates in reference cube
! misloc: location misfit abs(xyz - XYZ(uvw))
! is_inside: flag whether the target point locates inside the element
!

  real(kind=CUSTOM_REAL), intent(in) :: anchor_xyz(NDIM, NGNOD)
  real(kind=CUSTOM_REAL), intent(in) :: xyz(NDIM)

  real(kind=CUSTOM_REAL), intent(out) :: uvw(NDIM)
  real(kind=CUSTOM_REAL), intent(out) :: misloc 
  logical, intent(out) :: is_inside
  
  ! local variables
  integer :: iter_loop
  real(kind=CUSTOM_REAL), dimension(NDIM) :: xyzi ! iteratively improved xyz
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM) :: DuvwDxyz
  real(kind=CUSTOM_REAL), dimension(NDIM) :: dxyz, duvw

  real(kind=CUSTOM_REAL), parameter ::  zero = 0.0_CUSTOM_REAL & 
                                      , one = 1.0_CUSTOM_REAL &
                                      , minus_one = -1.0_CUSTOM_REAL

  ! initialize 
  uvw = zero
  is_inside = .true.

  ! iteratively update local coordinate uvw to approach the target xyz
  do iter_loop = 1, NUM_ITER

    ! predicted xyzi and Jacobian for the current uvw
    call cube2xyz(anchor_xyz, uvw, xyzi, DuvwDxyz)

    ! compute difference
    dxyz = xyz - xyzi

    ! compute increments
    duvw = matmul(DuvwDxyz, dxyz)

    ! update values
    uvw = uvw + duvw

    ! limit inside the cube
    if (any(uvw < minus_one .or. uvw > one)) then 
      where (uvw < minus_one) uvw = minus_one
      where (uvw > one) uvw = one
      ! flag inside in the last iteration step
      if (iter_loop == NUM_ITER) is_inside = .false.
    end if

  end do ! do iter_loop = 1,NUM_ITER
  
  ! calculate the predicted position 
  call cube2xyz(anchor_xyz, uvw, xyzi, DuvwDxyz)

  ! residual distance from the target point
  misloc = sqrt(sum((xyz-xyzi)**2))

end subroutine xyz2cube_bounded


!///////////////////////////////////////////////////////////////////////////////
subroutine cube2xyz(anchor_xyz, uvw, xyz, DuvwDxyz)
!-map from local coordinate (uvw) to physical position (xyz)
! the shape the element is defined by the anchor points (anchor_xyz)
!
!-input
! anchor_xyz(3): anchor points
! uvw(3): local coordinate 
!
!-output
! xyz(3): map uvw to physical space
! DuvwDxyz(3,3): jacobian matrix

  real(kind=CUSTOM_REAL), intent(in) :: anchor_xyz(NDIM, NGNOD)
  real(kind=CUSTOM_REAL), intent(in) :: uvw(NDIM)

  real(kind=CUSTOM_REAL), intent(out) :: xyz(NDIM)
  real(kind=CUSTOM_REAL), intent(out) :: DuvwDxyz(NDIM, NDIM)

  ! local variables
  real(kind=CUSTOM_REAL), dimension(NDIM) :: lag1, lag2, lag3 
  real(kind=CUSTOM_REAL), dimension(NDIM) :: lag1p, lag2p, lag3p
  real(kind=CUSTOM_REAL), dimension(NGNOD) :: shape3D 
  real(kind=CUSTOM_REAL), dimension(NGNOD, NDIM) :: dershape3D
  real(kind=CUSTOM_REAL) :: jacobian
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM) :: DxyzDuvw

  if (NGNOD /= 27) stop "ERROR: elements should have 27 control nodes"

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
  xyz = matmul(anchor_xyz, shape3D)
  DxyzDuvw = matmul(anchor_xyz, dershape3D)

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
  if (jacobian<=0.0) stop '3D Jacobian undefined'

  ! inverse matrix: Duvw/Dxyz = inv(Dxyz/Duvw) = adj(DxyzDuvw)/det(DxyzDuvw)  
  DuvwDxyz = DuvwDxyz/jacobian

end subroutine cube2xyz


!///////////////////////////////////////////////////////////////////////////////
subroutine hex_nodes(iaddx, iaddy, iaddz)
! define the topology of the hexahedral elements
! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif

! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

  if (NGNOD /= 27) stop 'elements should have 27 control nodes'

! corner nodes
  iaddx(1) = 0; iaddy(1) = 0; iaddz(1) = 0
  iaddx(2) = 2; iaddy(2) = 0; iaddz(2) = 0
  iaddx(3) = 2; iaddy(3) = 2; iaddz(3) = 0
  iaddx(4) = 0; iaddy(4) = 2; iaddz(4) = 0
  iaddx(5) = 0; iaddy(5) = 0; iaddz(5) = 2
  iaddx(6) = 2; iaddy(6) = 0; iaddz(6) = 2
  iaddx(7) = 2; iaddy(7) = 2; iaddz(7) = 2
  iaddx(8) = 0; iaddy(8) = 2; iaddz(8) = 2

! midside nodes (nodes located in the middle of an edge)
  iaddx(9) = 1; iaddy(9) = 0; iaddz(9) = 0
  iaddx(10) = 2; iaddy(10) = 1; iaddz(10) = 0
  iaddx(11) = 1; iaddy(11) = 2; iaddz(11) = 0
  iaddx(12) = 0; iaddy(12) = 1; iaddz(12) = 0
  iaddx(13) = 0; iaddy(13) = 0; iaddz(13) = 1
  iaddx(14) = 2; iaddy(14) = 0; iaddz(14) = 1
  iaddx(15) = 2; iaddy(15) = 2; iaddz(15) = 1
  iaddx(16) = 0; iaddy(16) = 2; iaddz(16) = 1
  iaddx(17) = 1; iaddy(17) = 0; iaddz(17) = 2
  iaddx(18) = 2; iaddy(18) = 1; iaddz(18) = 2
  iaddx(19) = 1; iaddy(19) = 2; iaddz(19) = 2
  iaddx(20) = 0; iaddy(20) = 1; iaddz(20) = 2

! side center nodes (nodes located in the middle of a face)
  iaddx(21) = 1; iaddy(21) = 1; iaddz(21) = 0
  iaddx(22) = 1; iaddy(22) = 0; iaddz(22) = 1
  iaddx(23) = 2; iaddy(23) = 1; iaddz(23) = 1
  iaddx(24) = 1; iaddy(24) = 2; iaddz(24) = 1
  iaddx(25) = 0; iaddy(25) = 1; iaddz(25) = 1
  iaddx(26) = 1; iaddy(26) = 1; iaddz(26) = 2

! center node (barycenter of the eight corners)
  iaddx(27) = 1; iaddy(27) = 1; iaddz(27) = 1

end subroutine hex_nodes


!///////////////////////////////////////////////////////////////////////////////
subroutine heap_sort(N,RA,IX)
!- Index an array RA of length N in ascending order by the Heapsort method
!
!-INPUTS:                                           
!      N     size of table RA                       
!      RA    table to be sorted                     
!-OUTPUT:                                           
!      IX    index table sorted in ascending order  
!                                                   
!-NOTE: The Heapsort method is a N Log2 N routine,  
!       and can be used for very large arrays.      

  integer, intent(in) :: N
  real(kind=CUSTOM_REAL), intent(in) :: RA(N)
  integer, intent(out) :: IX(N)

  ! local variables
  integer :: I, J, L, IR, IIX

  L=N/2+1
  IR=N
  do I=1,N
    IX(I) = I
  end do

  if (N<2) return 

  do

    if(L>1)then
      L=L-1
      IIX=IX(L)
    else
      IIX=IX(IR)
      IX(IR)=IX(1)
      IR=IR-1
      if(IR.eq.1)then
        IX(1)=IIX
        return
      end if
    end if
  
    I=L
    J=L+L
  
    do while (J <= IR)
      if (J < IR)then
        if(RA(IX(J)) < RA(IX(J+1)))  J=J+1
      end if
      if (RA(IIX) < RA(IX(J)))then
        IX(I)=IX(J)
        I=J; J=J+J
      else
        J=IR+1
      end if
    end do
  
    IX(I)=IIX
  end do

end subroutine heap_sort


end module sem_mesh
