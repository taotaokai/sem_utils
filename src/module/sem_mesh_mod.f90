module sem_mesh

!----------------------------------
! specification part
!----------------------------------

  use sem_constants, only: dp, CUSTOM_REAL, MAX_STRING_LEN, IIN, IOUT
  use sem_constants, only: GAUSSALPHA, GAUSSBETA
  use sem_constants, only: NGLLX, NGLLY, NGLLZ, NGNOD
! use sem_constants, only: ANGULAR_WIDTH_XI_IN_DEGREES_VAL
! use sem_constants, only: DEGREES_TO_RADIANS, R_UNIT_SPHERE
! use sem_constants, only: NEX_XI_VAL
! use sem_constants, only: NUM_ITER
  use sem_constants, only: R_EARTH_KM
  use sem_constants, only: REGIONAL_MOHO_MESH
  use sem_constants, only: IFLAG_CRUST, IFLAG_670_220, IFLAG_DUMMY

  use sem_io, only: sem_io_open_file_for_read

  implicit none

  private

  ! mesh data type
  type :: sem_mesh_data
    integer :: nspec, nglob
    real(dp), dimension(:,:), allocatable :: xyz_glob
    integer, dimension(:,:,:,:), allocatable :: ibool
    integer, dimension(:), allocatable :: idoubling
    !logical, dimension(:), allocatable :: ispec_is_tiso
    !real(dp), dimension(:,:,:,:), allocatable :: &
    !  xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  end type sem_mesh_data

  ! location within sem mesh data type 
  type :: sem_mesh_location
    integer :: stat
    integer :: eid
    real(dp) :: uvw(3)
    real(dp) :: misloc
    real(dp) :: lagrange(NGLLX, NGLLY, NGLLZ)
  end type sem_mesh_location

  ! public data types
  public :: sem_mesh_data
  public :: sem_mesh_location

  ! public operations
  public :: sem_mesh_init
  public :: sem_mesh_read
  !public :: sem_mesh_get_gll_volume
  !public :: sem_mesh_locate_xyz
  public :: sem_mesh_locate_kdtree2

!----------------------------------
! implementation part
!----------------------------------
contains

!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_init(mesh_data, nspec, nglob)
! allocate struct of mesh

  type(sem_mesh_data), intent(inout) :: mesh_data

  integer, intent(in) :: nspec, nglob

  if (allocated(mesh_data%xyz_glob)) then

    deallocate(mesh_data%xyz_glob, &
               mesh_data%ibool, &
               mesh_data%idoubling) !, &
               !mesh_data%ispec_is_tiso, &
               !mesh_data%xix, &
               !mesh_data%xiy, &
               !mesh_data%xiz, &
               !mesh_data%etax, &
               !mesh_data%etay, &
               !mesh_data%etaz, &
               !mesh_data%gammax, &
               !mesh_data%gammay, &
               !mesh_data%gammaz)
  end if

  allocate(mesh_data%xyz_glob(3, nglob), &
           mesh_data%ibool(NGLLX, NGLLY, NGLLZ, nspec), &
           mesh_data%idoubling(nspec)) !, &
           !mesh_data%ispec_is_tiso(nspec), &
           !   mesh_data%xix(NGLLX, NGLLY, NGLLZ, nspec), &
           !   mesh_data%xiy(NGLLX, NGLLY, NGLLZ, nspec), &
           !   mesh_data%xiz(NGLLX, NGLLY, NGLLZ, nspec), &
           !  mesh_data%etax(NGLLX, NGLLY, NGLLZ, nspec), &
           !  mesh_data%etay(NGLLX, NGLLY, NGLLZ, nspec), &
           !  mesh_data%etaz(NGLLX, NGLLY, NGLLZ, nspec), &
           !mesh_data%gammax(NGLLX, NGLLY, NGLLZ, nspec), &
           !mesh_data%gammay(NGLLX, NGLLY, NGLLZ, nspec), &
           !mesh_data%gammaz(NGLLX, NGLLY, NGLLZ, nspec))

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_read(mesh_data, basedir, iproc, iregion)
!-read mesh data (solver_data.bin)
!
!-inputs:
! (type mesh) mesh_data: structure mesh
! (string) basedir,iproc,iregion: specify the mesh file
! (integer) ngllxyz(3) = [NGLLX, NGLLY, NGLLZ]
!   gll points dimensions of the reference cube

  type(sem_mesh_data), intent(inout) :: mesh_data

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  integer :: nspec, nglob
  integer :: num, id, ispec, iglob
  real(dp) :: xyz_center(3), depth
  real(CUSTOM_REAL), allocatable :: dummy(:)

  call sem_io_open_file_for_read(IIN,basedir,iproc,iregion,'solver_data')

  ! get mesh dimensions
  read(IIN) nspec 
  read(IIN) nglob

  ! allocate memory for mesh_data
  call sem_mesh_init(mesh_data, nspec, nglob)

  mesh_data%nspec = nspec
  mesh_data%nglob = nglob

  allocate(dummy(nglob))

  read(IIN) dummy
  mesh_data%xyz_glob(1,:) = real(dummy, kind=dp)

  read(IIN) dummy
  mesh_data%xyz_glob(2,:) = real(dummy, kind=dp)

  read(IIN) dummy
  mesh_data%xyz_glob(3,:) = real(dummy, kind=dp)

  read(IIN) mesh_data%ibool
  read(IIN) mesh_data%idoubling
  !read(IIN) mesh_data%ispec_is_tiso

  !read(IIN) mesh_data%xix
  !read(IIN) mesh_data%xiy
  !read(IIN) mesh_data%xiz
  !read(IIN) mesh_data%etax
  !read(IIN) mesh_data%etay
  !read(IIN) mesh_data%etaz
  !read(IIN) mesh_data%gammax
  !read(IIN) mesh_data%gammay
  !read(IIN) mesh_data%gammaz

  close(IIN)

  ! separate crustal mesh layers for REGIONAL_MOHO_MESH 
  ! 3-layer crust: 10(third layer), 11, 12(shallowest layer)
  if (REGIONAL_MOHO_MESH) then
    num = 0
    do ispec = 1, nspec
      if (mesh_data%idoubling(ispec) == IFLAG_CRUST) then
        id = num - (num/3)*3
        mesh_data%idoubling(ispec) = 10 * IFLAG_CRUST + id
        num = num + 1
      endif
    enddo
  endif

  ! separate mesh layers across 410-km
  ! 40: above 410, 41: below 410
  do ispec = 1, nspec

    if (mesh_data%idoubling(ispec) == IFLAG_670_220) then

      iglob = mesh_data%ibool(NGLLX/2, NGLLY/2, NGLLZ/2, ispec)
      ! element center coordinate
      xyz_center = mesh_data%xyz_glob(:,iglob)
      depth = (1.0 - sqrt(sum(mesh_data%xyz_glob(:,iglob)**2))) * R_EARTH_KM

      if (depth < 405.0_CUSTOM_REAL) then
        mesh_data%idoubling(ispec) = 10 * IFLAG_670_220
      else ! below 410-km
        mesh_data%idoubling(ispec) = 10 * IFLAG_670_220 + 1
      endif

    endif

  enddo

end subroutine


!///////////////////////////////////////////////////////////////////////////////
!subroutine sem_mesh_get_gll_volume(mesh_data, gll_volume)
!
!! compute volume corresponding to each gll point in each spectral element
!! gll_volume = jacobian * wgll_cube
!! jacobian = det(d(x,y,z)/d(xi,eta,gamma)), volumetric deformation ratio from cube to actual element 
!
!  type (mesh), intent(in) :: mesh_data
!  real(kind=CUSTOM_REAL), intent(out) :: gll_volume(NGLLX,NGLLY,NGLLZ,NSPEC)
!
!  integer :: i,j,k,ispec
!  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
!  real(kind=CUSTOM_REAL) :: jacobianl
!
!  ! calculate gll quad weight in the reference cube
!  !call compute_wgll_cube()
!
!  ! calculate volume integration weight on gll
!  do ispec = 1, NSPEC
!    do k = 1, NGLLZ
!      do j = 1, NGLLY
!        do i = 1, NGLLX
!
!          ! gets derivatives of ux, uy and uz with respect to x, y and z
!                xixl = mesh_data%xix(i,j,k,ispec)
!                xiyl = mesh_data%xiy(i,j,k,ispec)
!                xizl = mesh_data%xiz(i,j,k,ispec)
!              etaxl = mesh_data%etax(i,j,k,ispec)
!              etayl = mesh_data%etay(i,j,k,ispec)
!              etazl = mesh_data%etaz(i,j,k,ispec)
!          gammaxl = mesh_data%gammax(i,j,k,ispec)
!          gammayl = mesh_data%gammay(i,j,k,ispec)
!          gammazl = mesh_data%gammaz(i,j,k,ispec)
!
!          ! jacobian: det( d(x,y,z)/d(xi,eta,gamma))
!          jacobianl = 1.0_CUSTOM_REAL / &
!            ( xixl*(etayl*gammazl-etazl*gammayl)  &
!             -xiyl*(etaxl*gammazl-etazl*gammaxl)  &
!             +xizl*(etaxl*gammayl-etayl*gammaxl))
!
!          ! gll volume: jacobian * wgll_cube 
!          gll_volume(i,j,k,ispec) = jacobianl * wgll_cube(i,j,k)
!
!        enddo
!      enddo
!    enddo
!  enddo
!  
!end subroutine


!///////////////////////////////////////////////////////////////////////////////
!subroutine compute_wgll_cube()
!! compute the integration weight for gll quadrature in a reference cube
!! i.e. volume corresponding to gll point in the reference cube
!
!  integer :: i, j, k
!
!  call zwgljd(xigll, wxgll, NGLLX, GAUSSALPHA, GAUSSBETA)
!  call zwgljd(yigll, wygll, NGLLY, GAUSSALPHA, GAUSSBETA)
!  call zwgljd(zigll, wzgll, NGLLZ, GAUSSALPHA, GAUSSBETA)
!
!  do k=1, NGLLZ
!    do j=1, NGLLY
!      do i=1, NGLLX
!        wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
!      enddo
!    enddo
!  enddo
!
!end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_mesh_locate_kdtree2( &
  mesh_data, npoint, xyz, idoubling, nnearest, max_misloc, &
  location_result)
!-locate xyz(3,npoint) inside a given SEM mesh
!
!-inputs:
! type(sem_mesh_data) mesh_data
! (real) xyz(3, npoint): target points to locate
! (int) idoubling(npoint): layer id of targt points 
!   not specified if = IFLAG_DUMMY, see sem_mesh_read() for more
! (int) nnearest: number of nearest global points in SEM mesh for kdtree search
! (real) max_misloc: mislocation threshold
!
!-output
! type(sem_mesh_location) location_result(npoint): location results

  use kdtree2_module

  type(sem_mesh_data), intent(in) :: mesh_data
  integer, intent(in) :: npoint
  real(dp), intent(in) :: xyz(3, npoint)
  integer, intent(in) :: idoubling(npoint)
  integer, intent(in) :: nnearest
  real(dp), intent(in) :: max_misloc

  type(sem_mesh_location), intent(out) :: location_result(npoint) 

  !-- local varaibles
  integer :: nspec
  integer :: i, ipoint, ispec, iglob

  !-- anchor points
  integer :: ia
  integer, dimension(NGNOD) :: iax, iay, iaz 
  real(dp), allocatable :: xyz_anchor(:,:,:)

  !-- kdtree2
  type(kdtree2), pointer :: tree
  type(kdtree2_result), allocatable :: search_result(:)

  !-- find element 
  integer :: inn, isel, nsel
  real(dp) :: two_max_misloc
  integer, allocatable :: RANGE_1_NSPEC(:), eid_sel(:)
  logical, allocatable :: flag_ibool(:,:,:,:), flag_nspec(:)

  !-- locate inside element
  integer, parameter :: niter_xyz2cube = 4
  real(dp) :: xyz1(3), uvw1(3), misloc1
  logical :: flag_inside

  !-- GLL colocation points and lagrange interpolation weights
  real(dp), dimension(NGLLX) :: xigll, wxgll
  real(dp), dimension(NGLLY) :: yigll, wygll
  real(dp), dimension(NGLLZ) :: zigll, wzgll
  real(dp) :: hlagx(NGLLX), hlagy(NGLLY), hlagz(NGLLZ)
  integer :: igllx, iglly, igllz

  !===== mesh dimension variable

  nspec = mesh_data%nspec

  !===== coordinates of GLL points

  call zwgljd(xigll, wxgll, NGLLX, GAUSSALPHA, GAUSSBETA)
  call zwgljd(yigll, wygll, NGLLY, GAUSSALPHA, GAUSSBETA)
  call zwgljd(zigll, wzgll, NGLLZ, GAUSSALPHA, GAUSSBETA)

  !===== anchor points of mesh elements 

  ! get index of anchor points in the GLL element
  call anchor_point_index(iax, iay, iaz)

  ! get anchor and center points of each GLL element 
  allocate(xyz_anchor(3, NGNOD, nspec))
  do ispec = 1, nspec
    do ia = 1, NGNOD
      iglob = mesh_data%ibool(iax(ia), iay(ia), iaz(ia), ispec)
      xyz_anchor(:, ia, ispec) = mesh_data%xyz_glob(:, iglob)
    enddo
    ! the last anchor point is the element center
    !xyz_center(:,ispec) = mesh_data%xyz(:,iglob)
  enddo

  !===== build kdtree

  tree => kdtree2_create(mesh_data%xyz_glob, sort=.true., rearrange=.true.)

  !===== locate each point

  !-- intialize variables and arrays
  two_max_misloc = 2.0 * max_misloc

  allocate(search_result(nnearest))
  allocate(RANGE_1_NSPEC(nspec), eid_sel(nspec), flag_nspec(nspec))
  allocate(flag_ibool(NGLLX, NGLLY, NGLLZ, nspec))

  RANGE_1_NSPEC = [ (i, i=1,nspec) ]

  ! initialize location results
  location_result(:)%stat = -1
  location_result(:)%eid = -1
  location_result(:)%misloc = huge(1.0_dp)

  loop_points: do ipoint = 1, npoint

    xyz1 = xyz(:,ipoint)
    
    !-- get the n nearest points in the mesh
    call kdtree2_n_nearest(tp=tree, qv=xyz1, nn=nnearest, &
      results=search_result)

    !-- find the elements associated with the n nearest points
    flag_ibool = .false.
    do inn = 1, nnearest
      if (search_result(inn)%dis > two_max_misloc) then
        cycle
      endif
      iglob = search_result(inn)%idx
      flag_ibool = flag_ibool .or. (mesh_data%ibool==iglob)
    enddo

    flag_nspec = .false.
    do ispec = 1, nspec
      if (any(flag_ibool(:,:,:,ispec))) then
        flag_nspec(ispec) = .true.
      endif
    enddo

    !-- test each potential element 
    nsel = count(flag_nspec)
    eid_sel = pack(RANGE_1_NSPEC, mask=flag_nspec)

    do isel = 1, nsel

      ispec = eid_sel(isel)

      ! test layer id (if used)
      if (idoubling(ipoint) /= IFLAG_DUMMY .and. &
          idoubling(ipoint) /= mesh_data%idoubling(ispec)) then
        cycle
      endif

      ! locate point to this element
      call xyz2cube_bounded(xyz_anchor(:,:,ispec), xyz1, niter_xyz2cube, &
        uvw1, misloc1, flag_inside)

      if (flag_inside) then ! record this element and exit looping the rest elements
        location_result(ipoint)%stat = 1
        location_result(ipoint)%eid = ispec
        location_result(ipoint)%misloc = misloc1
        location_result(ipoint)%uvw = uvw1
        exit

      else ! accept this element if misloc smaller than 
           ! max_misloc and previous misloc

        if (misloc1 < max_misloc .and. &
            misloc1 < location_result(ipoint)%misloc) then

          location_result(ipoint)%stat = 0
          location_result(ipoint)%eid = ispec
          location_result(ipoint)%misloc = misloc1
          location_result(ipoint)%uvw = uvw1
        endif

      endif ! flag_inside

    enddo ! isel

    ! set interpolation weights on GLL points if located
    if (location_result(ipoint)%stat /= -1) then

      call lagrange_poly(location_result(ipoint)%uvw(1), NGLLX, xigll, hlagx)
      call lagrange_poly(location_result(ipoint)%uvw(2), NGLLY, yigll, hlagy)
      call lagrange_poly(location_result(ipoint)%uvw(3), NGLLZ, zigll, hlagz)

      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            location_result(ipoint)%lagrange(igllx,iglly,igllz) = &
              hlagx(igllx) * hlagy(iglly) * hlagz(igllz)
          enddo
        enddo
      enddo 

    endif

  enddo loop_points 

end subroutine


!///////////////////////////////////////////////////////////////////////////////
!subroutine sem_mesh_locate_xyz( &
!  mesh_data, npoint, xyz, idoubling, num_typical_size, &
!  uvw, hlagrange, misloc, eid, locstat)
!!-locate xyz(3,npoint) inside a given SEM mesh
!!
!!-inputs:
!! type(mesh) mesh_data
!! (real) xyz(3, npoint): target points to locate
!! (int) idoubling(npoint): layer id of targt points 
!!   not specified if = IFLAG_DUMMY, see sem_mesh_read() for more
!! (int) num_typical_size: number of typical element size for locating points
!!
!!-output
!! uvw(3, npoint): local coordinates inside the located element
!! hlagrange(NGLLX, NGLLY, NGLLZ, npoint): GLL interpolation weight at each GLL
!!   point inside the element
!! misloc(npoint): distance between the target and predicted position abs(xyz - XYZ(uvw))
!! eid(npoint): id of the element which contains the target point
!! locstat(npoint): status of location 
!!   (1: inside, 0: close to the element, -1: out)
!
!  type (mesh), intent(in) :: mesh_data
!  integer, intent(in) :: npoint
!  real(dp), intent(in) :: xyz(3, npoint)
!  integer, intent(in) :: idoubling(npoint)
!  integer, intent(in) :: num_typical_size
!
!  real(dp), intent(out) :: uvw(3, npoint)
!  real(dp), intent(out) :: hlagrange(NGLLX, NGLLY, NGLLZ, npoint)
!  real(dp), intent(out) :: misloc(npoint)
!  integer, intent(out) :: eid(npoint)
!  integer, intent(out) :: locstat(npoint)
!
!  ! local varaibles
!  integer :: nspec, nglob
!  integer :: i, ia, nsel, ispec, iglob, ipoint, isel
!  integer :: RANGE_1_NSPEC(nspec), eid_sel(nspec), ind_sort(nspec)
!  real(dp) :: typical_size2, typical_size, dist_sq(nspec)
!  real(dp) :: max_x, min_x, max_y, min_y, max_z, min_z
!  real(dp) :: xyz1(3), uvw1(3), misloc1
!  logical :: idx_sel(nspec), is_inside
!  integer :: igllx, iglly, igllz
!  real(dp) :: hlagx(NGLLX), hlagy(NGLLY), hlagz(NGLLZ)
!  !-- index of anchor points inside the spectral element
!  integer, dimension(NGNOD) :: iax, iay, iaz 
!  !-- center and anchor points inside each spectral element
!  real(dp) :: center_xyz(3, NSPEC)
!  real(dp) :: anchor_xyz(3, NGNOD, NSPEC)
!
!  !-- GLL points
!  real(dp), dimension(NGLLX) :: xigll, wxgll
!  real(dp), dimension(NGLLY) :: yigll, wygll
!  real(dp), dimension(NGLLZ) :: zigll, wzgll
!
!  !==== initialization
!
!  ! get index of anchor points in the GLL element
!  call anchor_point_index(iax,iay,iaz)
!
!  ! get anchor and center points of each GLL element 
!  do ispec = 1, NSPEC
!    do ia = 1, NGNOD
!      iglob = mesh_data%ibool(iax(ia),iay(ia),iaz(ia),ispec)
!      anchor_xyz(:,ia,ispec) = mesh_data%xyz(:,iglob)
!    enddo
!    ! the last anchor point is the element center
!    center_xyz(:,ispec) = mesh_data%xyz(:,iglob)
!  enddo
!
!  ! initialize some parameters for locating
!  RANGE_1_NSPEC = [ (i, i=1,nspec) ]
!  uvw = huge(1.0_dp)
!  hlagrange = huge(1.0_dp)
!  misloc = huge(1.0_dp)
!  locstat = -1
!  eid = -1
!  ! search range: typical element size at surface
!  typical_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
!                     ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
!                 * DEGREES_TO_RADIANS * R_UNIT_SPHERE
!  typical_size2 = typical_size * num_typical_size
!  ! minimum/maximum x,y,z coordinates in the mesh
!  max_x = maxval(center_xyz(1,:)) + typical_size2
!  min_x = minval(center_xyz(1,:)) - typical_size2
!  max_y = maxval(center_xyz(2,:)) + typical_size2
!  min_y = minval(center_xyz(2,:)) - typical_size2
!  max_z = maxval(center_xyz(3,:)) + typical_size2
!  min_z = minval(center_xyz(3,:)) - typical_size2
!
!  ! coordinates of GLL points
!  call zwgljd(xigll, wxgll, NGLLX, GAUSSALPHA, GAUSSBETA)
!  call zwgljd(yigll, wygll, NGLLY, GAUSSALPHA, GAUSSBETA)
!  call zwgljd(zigll, wzgll, NGLLZ, GAUSSALPHA, GAUSSBETA)
!
!  !===== locate each point
!
!  do ipoint = 1, npoint
!
!    ! target point
!    xyz1 = xyz(:,ipoint)
!
!    ! skip points based on coordinate ranges
!    if (     xyz1(1)<min_x .or. xyz1(1)>max_x &
!        .or. xyz1(2)<min_y .or. xyz1(2)>max_y &
!        .or. xyz1(3)<min_z .or. xyz1(3)>max_z) then
!      cycle
!    end if
!
!    ! select elements of the same layer id (idoubling)
!    idx_sel = .true.
!    if (idoubling(ipoint) /= IFLAG_DUMMY) then
!      idx_sel = mesh_data%idoubling == idoubling(ipoint)
!      if (count(idx_sel) == 0) then 
!        cycle
!      endif
!    endif
!    ! also select elements close to the target point within a certain range
!    idx_sel = idx_sel .and. &
!              (center_xyz(1,:)>xyz1(1)-typical_size2) .and. &
!              (center_xyz(1,:)<xyz1(1)+typical_size2) .and. &
!              (center_xyz(2,:)>xyz1(2)-typical_size2) .and. &
!              (center_xyz(2,:)<xyz1(2)+typical_size2) .and. &
!              (center_xyz(3,:)>xyz1(3)-typical_size2) .and. &
!              (center_xyz(3,:)<xyz1(3)+typical_size2)
!    nsel = count(idx_sel)
!    if (nsel == 0) then 
!      cycle
!    endif
!
!    ! sort selected elements based on distaces to the target point
!    eid_sel(1:nsel) = pack(RANGE_1_NSPEC, idx_sel)
!    dist_sq(1:nsel) =  (center_xyz(1,eid_sel(1:nsel)) - xyz1(1))**2 & 
!                     + (center_xyz(2,eid_sel(1:nsel)) - xyz1(2))**2 &
!                     + (center_xyz(3,eid_sel(1:nsel)) - xyz1(3))**2
!
!    call heap_sort(nsel, dist_sq, ind_sort)
!
!    ! loop all selected elements to find the one containing the target xyz 
!    do isel = 1, nsel
!
!      ispec = eid_sel(ind_sort(isel))
!
!      call xyz2cube_bounded(anchor_xyz(:,:,ispec),xyz1,uvw1,misloc1,is_inside)
!
!      if (is_inside) then ! record this element and exit looping the rest elements
!        locstat(ipoint) = 1
!        eid(ipoint) = ispec
!        misloc(ipoint) = misloc1
!        uvw(:,ipoint) = uvw1
!        exit
!      else ! record this element if smaller misloc
!        if (misloc1 < misloc(ipoint)) then
!          locstat(ipoint) = 0
!          eid(ipoint) = ispec
!          misloc(ipoint) = misloc1
!          uvw(:,ipoint) = uvw1
!        endif
!      endif !inside
!
!    enddo ! isel
!
!    ! set interpolation weights on GLL points if located
!    if (locstat(ipoint) /= -1) then
!
!      call lagrange_poly(uvw(1,ipoint), NGLLX, xigll, hlagx)
!      call lagrange_poly(uvw(2,ipoint), NGLLY, yigll, hlagy)
!      call lagrange_poly(uvw(3,ipoint), NGLLZ, zigll, hlagz)
!
!      do igllz = 1, NGLLZ
!        do iglly = 1, NGLLY
!          do igllx = 1, NGLLX
!            hlagrange(igllx,iglly,igllz,ipoint) &
!              = hlagx(igllx) * hlagy(iglly) * hlagz(igllz)
!          enddo
!        enddo
!      enddo 
!
!    endif
!
!  enddo ! ipoint
!
!end subroutine


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

  real(dp), intent(in) :: x
  integer, intent(in) :: ngll
  real(dp), intent(in) :: xgll(ngll)

  real(dp), intent(out) :: lagrange(ngll)

  ! local variables
  integer :: i
  integer, dimension(ngll) :: ind
  real(dp), dimension(ngll) :: xx, yy

  ! lagrange(ngll) 
  ind = (/(i, i=1,ngll)/)
  xx = x - xgll

  do i = 1, ngll
    yy = xgll(i) - xgll
    lagrange(i) = product(xx/yy, mask=(ind/=i))
  end do

end subroutine lagrange_poly

!///////////////////////////////////////////////////////////////////////////////
subroutine xyz2cube_bounded(xyz_anchor, xyz, niter, uvw, misloc, flag_inside )
!-mapping a given point in physical space (xyz) to the 
! reference cube (uvw), 
! and also flag whether the point is inside the cube
! if the point lies outside the element, calculate the bounded (xi,eta,gamma)
! inside or on the surface of the reference unit cube.
!
!-inputs:
! (real) xyz_anchor(3, NGNOD): anchor points of the element
! xyz: coordinates of the target point
! (int) niter: number of iterations to invert for uvw
!
!-outputs:
! uvw(3): local coordinates in reference cube
! misloc: location misfit abs(xyz - XYZ(uvw))
! flag_inside: flag whether the target point locates inside the element

  real(dp), intent(in) :: xyz_anchor(3, NGNOD)
  real(dp), intent(in) :: xyz(3)
  integer, intent(in) :: niter

  real(dp), intent(out) :: uvw(3)
  real(dp), intent(out) :: misloc 
  logical, intent(out) :: flag_inside
  
  ! local variables
  integer :: iter
  real(dp), dimension(3) :: xyzi ! iteratively improved xyz
  real(dp), dimension(3,3) :: DuvwDxyz
  real(dp), dimension(3) :: dxyz, duvw

  real(dp), parameter ::  zero = 0.0_dp, & 
                          one = 1.0_dp, &
                          minus_one = -1.0_dp

  ! initialize 
  uvw = zero
  flag_inside = .true.

  ! iteratively update local coordinate uvw to approach the target xyz
  do iter = 1, niter

    ! predicted xyzi and Jacobian for the current uvw
    call cube2xyz(xyz_anchor, uvw, xyzi, DuvwDxyz)

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
      ! set is_inside to false based on the last iteration
      if (iter == niter) then
        flag_inside = .false.
      endif
    endif

  enddo ! do iter_loop = 1,NUM_ITER
  
  ! calculate the predicted position 
  call cube2xyz(xyz_anchor, uvw, xyzi, DuvwDxyz)

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

  real(dp), intent(in) :: anchor_xyz(3, NGNOD)
  real(dp), intent(in) :: uvw(3)

  real(dp), intent(out) :: xyz(3)
  real(dp), intent(out) :: DuvwDxyz(3, 3)

  ! local variables
  real(dp), dimension(3) :: lag1, lag2, lag3 
  real(dp), dimension(3) :: lag1p, lag2p, lag3p
  real(dp), dimension(NGNOD) :: shape3D 
  real(dp), dimension(NGNOD, 3) :: dershape3D
  real(dp) :: jacobian
  real(dp), dimension(3,3) :: DxyzDuvw

  if (NGNOD /= 27) then
    stop "[ERROR] cube2xyz: elements should have 27 control nodes"
  endif

  ! lagrange polynomials of order 3 on [-1,1], with collocation points: -1,0,1 
  lag1 = uvw * (uvw - 1.0_dp) / 2.0_dp
  lag2 = 1.0_dp - uvw**2
  lag3 = uvw * (uvw + 1.0_dp) / 2.0_dp
  
  ! derivative of lagrange polynomials
  lag1p = uvw - 0.5_dp
  lag2p = -2.0_dp * uvw
  lag3p = uvw + 0.5_dp
  
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
  if (jacobian<=0.0_dp) then
    print *, "[ERROR] cube2xyz: 3D Jacobian undefined jacobian=", jacobian
    stop
  endif

  ! inverse matrix: Duvw/Dxyz = inv(Dxyz/Duvw) = adj(DxyzDuvw)/det(DxyzDuvw)  
  DuvwDxyz = DuvwDxyz / jacobian

end subroutine cube2xyz


!///////////////////////////////////////////////////////////////////////////////
subroutine hex_nodes(iaddx, iaddy, iaddz)
! index of gll points that define the shape of element
! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif

  integer, dimension(NGNOD), intent(out) :: iaddx,iaddy,iaddz

  if (NGNOD /= 27) then
    stop '[ERROR] hex_nodes: elements should have 27 control nodes'
  endif

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
subroutine anchor_point_index(iax,iay,iaz)
!-get the index of anchor points in the GLL points of spectral element
! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif
! currently anchor points are defined on 3x3x3 hex
!
!-output: iax, iay, iaz

  ! index of anchor points
  integer, dimension(NGNOD), intent(out) :: iax,iay,iaz

  ! topology of the control points of the surface element
  integer :: ia
  integer, dimension(NGNOD) :: iaddx, iaddy, iaddz

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddz)

  do ia = 1, NGNOD

    if (iaddx(ia) == 0) then
      iax(ia) = 1
    else if (iaddx(ia) == 1) then
      iax(ia) = (NGLLX+1)/2
    else if (iaddx(ia) == 2) then
      iax(ia) = NGLLX
    else
      stop 'incorrect value of iaddx'
    endif

    if (iaddy(ia) == 0) then
      iay(ia) = 1
    else if (iaddy(ia) == 1) then
      iay(ia) = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay(ia) = NGLLY
    else
      stop 'incorrect value of iaddy'
    endif

    if (iaddz(ia) == 0) then
      iaz(ia) = 1
    else if (iaddz(ia) == 1) then
      iaz(ia) = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
      iaz(ia) = NGLLZ
    else
      stop 'incorrect value of iaddz'
    endif

  end do ! do ia = 1,NGNOD

end subroutine anchor_point_index


end module sem_mesh
