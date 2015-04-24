! interpolate the GLL model at a given location

subroutine read_mesh ( &
  datdir, iproc, NSPEC_MAX, NGLOB_MAX, &
  nspec, nglob, xyz, xyza, xyzc, ibool, idoubling, model )
!input:
!  datdir, iproc: database directory, slice id, region id
!output:
!  nspec, nglob: dimensions of mesh
!  xyz(nglob): coordinates of GLL points
!  xyza(3,na,nspec): anchor points defining element shape function 
!  xyzc(3,nspec): element center coordinates
!  ibool(NGLLX,NGLLY,NGLLZ,nspec): map GLL point in each element to global index (1:nglob) 
!  idoubling(nspec): mesh layer ID for each element
!  model(NGLLX,NGLLY,NGLLZ,nspec,nparam): model values on gll point in each element

  use constants, only: MAX_STRING_LEN,NDIM,NGLLX,NGLLY,NGLLZ,NGNOD,SIZE_REAL,IIN
  use user_par, only: RK,NPARS,PARNAME,IREG

  implicit none

  !==========================================================
  ! subroutine arguments 
  !==========================================================
  ! inputs
  character(len=MAX_STRING_LEN) :: datdir
  integer :: iproc, NSPEC_MAX, NGLOB_MAX

  ! outputs: mesh data
  integer :: nspec, nglob
  real(RK), dimension(NDIM,NGLOB_MAX) :: xyz
  real(RK), dimension(NDIM,NGNOD,NSPEC_MAX) :: xyza
  real(RK), dimension(NDIM,NSPEC_MAX) :: xyzc
  real(RK), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX,NPARS) :: model
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: ibool
  integer, dimension(NSPEC_MAX) :: idoubling

  !==========================================================
  ! local variables
  !==========================================================
  character(len=MAX_STRING_LEN) :: basename, fn_topo, fn_model
  integer :: ia,ispec,iglob,ipar
  integer,dimension(NGNOD) :: iax, iay, iaz ! index of anchor points

  !integer :: id, num, ier
  integer :: ier
  !real(8) :: depth

  real(SIZE_REAL), dimension(NDIM,NGLOB_MAX) :: xyz_in
  real(SIZE_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX,NPARS) :: model_in

  !==========================================================
  ! read topology data
  !==========================================================
  ! file basename
  !print *,'datdir,iproc,ireg:',trim(datdir),iproc,IREG
  write(basename,'(a,a,i6.6,a,i1,a)') trim(datdir),'/proc',iproc,'_reg',ireg,'_'

  ! open topology file
  fn_topo = trim(basename)//'solver_data.bin'
  open(unit=IIN, file=trim(fn_topo), status='old', action='read', iostat=ier, form='unformatted')
  if (ier /= 0) then
    print *,'file: ',trim(fn_topo)
    stop 'Error opening file'
  end if

  ! read mesh dimension
  read(IIN) nspec
  read(IIN) nglob

  !print *,'nspec,nglob',nspec,nglob
  if (nspec>NSPEC_MAX .or. nglob>NGLOB_MAX) stop 'read_mesh: NSPEC/NGLOB_MAX too small'

  ! read in mesh topology
  read(IIN) xyz_in(1,1:nglob)
  read(IIN) xyz_in(2,1:nglob)
  read(IIN) xyz_in(3,1:nglob)
  read(IIN) ibool(1:NGLLX,1:NGLLY,1:NGLLZ,1:nspec)
  read(IIN) idoubling(1:nspec)
  close(IIN)

  xyz(:,1:nglob) = DBLE(xyz_in(:,1:nglob))

  !==========================================================
  ! anchor/center points 
  !==========================================================
  ! index of anchor points within GLL points of a spectral element
  call anchor_point_index(iax,iay,iaz)

  ! get anchor/center grids 
  do ispec = 1,nspec
    do ia = 1,NGNOD
      iglob = ibool(iax(ia),iay(ia),iaz(ia),ispec)
      xyza(:,ia,ispec) = xyz(:,iglob)
    end do
    xyzc(:,ispec) = xyz(:,iglob) ! the last anchor point is the body center
  end do

  !==========================================================
  ! modify idoubling
  !==========================================================
! ! crust layering
! num = 0
! do ispec = 1,nspec
!   if (idoubling(ispec)==IFLAG_CRUST) then
!     id = num-(num/3)*3
!     idoubling(ispec) = 10*IFLAG_CRUST+NLAYER_CRUST-1-id
!     num = num+1
!   end if
! end do

! ! separate mesh layer across 410km
! do ispec = 1,nspec
!   if (idoubling(ispec)==IFLAG_670_220) then
!     depth = (1-sqrt(sum(xyzc(:,ispec)**2)))*R_EARTH_KM
!     if (depth < 405.d0) then
!       idoubling(ispec) = 10*IFLAG_670_220
!     else
!       idoubling(ispec) = 10*IFLAG_670_220+1
!     end if
!   end if
! end do

  !==========================================================
  ! read model data
  !==========================================================
  do ipar = 1,NPARS
    fn_model = trim(basename)//trim(PARNAME(ipar))//'.bin'

    open(unit=IIN, file=trim(fn_model), status='old', action='read', iostat=ier, form='unformatted')
    if (ier /= 0) then
      print *,'file: ',trim(fn_model)
      stop 'Error opening file'
    end if
    
    read(IIN,iostat=ier) model_in(:,:,:,1:nspec,ipar)
    if (ier /= 0) then
      print *,'file: ',trim(fn_model)
      stop 'Error reading model data'
    end if

    model(:,:,:,1:nspec,ipar) = DBLE(model_in(:,:,:,1:nspec,ipar)) 

    close(IIN)
  end do

end subroutine read_mesh
