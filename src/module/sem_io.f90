module sem_io

!----------------------------------
! specification part
!----------------------------------

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA, &
    IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE

  implicit none

private

  ! array dimension
  integer :: NSPEC, NGLOB
  logical :: is_dimension_set=.false.

  ! for calculation of volume correspoding to each gll point
  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll
  ! array with all the integration weights in the cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  logical :: is_wgll_cube_set = .false.

  ! mesh structure
  type, public :: sem_mesh
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyz
    integer, dimension(:,:,:,:), allocatable :: ibool
    integer, dimension(:), allocatable :: idoubling
    logical, dimension(:), allocatable :: ispec_is_tiso 
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
    ! volume corresponding to each gll point
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: volume_gll
    ! flag arrays allocated
    logical :: is_allocated = .false.
  end type sem_mesh

  !---- public
  ! dimensions 
  public :: NGLLX,NGLLY,NGLLZ,NSPEC,NGLOB

  ! subroutines
  public :: sem_set_dimension
  public :: sem_open_datafile_for_read
  public :: sem_open_datafile_for_write
  public :: sem_read_mesh
  public :: sem_set_volume_gll
  public :: sem_read_kernel_cijkl
  public :: sem_read_model
  public :: sem_read_model_n
  public :: sem_write_model
  public :: sem_write_model_n

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
subroutine sem_set_dimension(iregion)

  integer, intent(in) :: iregion

  select case (iregion)
    case (IREGION_CRUST_MANTLE)
      NSPEC = NSPEC_CRUST_MANTLE
      NGLOB = NGLOB_CRUST_MANTLE
!   case (IREGION_OUTER_CORE)
!     NSPEC = NSPEC_CRUST_MANTLE
!     NGLOB = NGLOB_CRUST_MANTLE
!   case (IREGION_INNER_CORE)
!     NSPEC = NSPEC_CRUST_MANTLE
!     NGLOB = NGLOB_CRUST_MANTLE
    case default
      write(*,*) 'sem_IO: unrecognized region code ( ', iregion, ' )'
      stop
  end select

  is_dimension_set = .true.

end subroutine

!//////////////////////////
subroutine sem_open_datafile_for_read(IIN,basedir,iproc,iregion,model_name)

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion, IIN

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(model_name)
  open(IIN,file=trim(filename),status='old', &
    form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'sem_IO: Error open file for read: ', trim(filename)
    stop
  endif

end subroutine

!//////////////////////////
subroutine sem_open_datafile_for_write(IOUT,basedir,iproc,iregion,model_name)

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion, IOUT

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(model_name)
  open(IOUT,file=trim(filename),status='unknown', &
    form='unformatted',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'sem_IO: Error open file for write: ', trim(filename)
    stop
  endif

end subroutine

!//////////////////////////
subroutine sem_init_mesh(mesh_data)

  type(sem_mesh), intent(inout) :: mesh_data

  if (.not. is_dimension_set) then
    print *, 'sem_IO: undetermined dimension, invoke sem_set_dimension() first!'
    stop
  end if

  if (mesh_data%is_allocated) then
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
               mesh_data%gammaz, &
               mesh_data%volume_gll)
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
               mesh_data%gammaz(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC), &
           mesh_data%volume_gll(1:NGLLX,1:NGLLY,1:NGLLZ,1:NSPEC))

  mesh_data%is_allocated = .true.

end subroutine

!//////////////////////////
subroutine sem_read_mesh(mesh_data, basedir, iproc, iregion)

  type(sem_mesh), intent(inout) :: mesh_data
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  integer :: ival

  if (.not. is_dimension_set) then
    print *, 'sem_IO: undetermined dimension, involk sem_set_dimension() first!'
    stop
  end if
  if (.not. mesh_data%is_allocated) then
    call sem_init_mesh(mesh_data)
  end if

  call sem_open_datafile_for_read(IIN,basedir,iproc,iregion,'solver_data')

  read(IIN) ival ! nspec
  if (ival /= NSPEC) then
    write(*,*) 'sem_IO: nspec=',ival,' in solver_data.bin, while NSPEC=',NSPEC
    stop
  end if

  read(IIN) ival ! nglob
  if (ival /= NGLOB) then
    write(*,*) 'sem_IO: nglob=',ival,' in solver_data.bin, while NGLOB=',NGLOB
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

!//////////////////////////
! compute the integration weight for gll quadrature in a reference cube
! i.e. volume corresponding to gll point in the reference cube
subroutine compute_wgll_cube()

  integer :: i,j,k

  ! GLL points
  wgll_cube = 0.0d0
  call _zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call _zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call _zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = real(wxgll(i)*wygll(j)*wzgll(k), kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  is_wgll_cube_set = .true.

end subroutine

!//////////////////////////
! compute volume corresponding to each gll point in each spectral element
! volume_gll = jacobian * wgll_cube
! jacobian = det(d(x,y,z)/d(xi,eta,gamma)), volumetric deformation ratio from cube to actual element 
subroutine sem_set_volume_gll(mesh_data)

  type(sem_mesh), intent(inout) :: mesh_data

  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: jacobianl

  ! calculate gll quad weight in the reference cube
  if (.not. is_wgll_cube_set) then
    call compute_wgll_cube()
  end if

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
          mesh_data%volume_gll(i,j,k,ispec) = jacobianl * wgll_cube(i,j,k)

        enddo
      enddo
    enddo
  enddo
  
end subroutine

!//////////////////////////
subroutine sem_read_model(model_array, basedir, iproc, iregion, model_name)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  integer :: ier

  if (.not. is_dimension_set) then
    print *, 'sem_IO: undetermined dimension, involk sem_set_dimension() first!'
    stop
  end if

  call sem_open_datafile_for_read(IIN,basedir,iproc,iregion,model_name)
  read(IIN, iostat=ier) model_array
  if (ier /= 0) stop "sem_read_model:read error"

  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_read_model_n(model_array, basedir, iproc, iregion, nmodel, model_names)

  integer, intent(in) :: nmodel
  real(kind=CUSTOM_REAL), dimension(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion

  integer :: i

  do i = 1, nmodel
    call sem_read_model(model_array(i,:,:,:,:), basedir, iproc, iregion, model_names(i))
  end do

end subroutine

!//////////////////////////
subroutine sem_read_kernel_cijkl(model_array, basedir, iproc, iregion)

  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  if (.not. is_dimension_set) then
    print *, 'sem_IO: undetermined dimension, involk set_dimension() first!'
    stop
  end if

  call sem_open_datafile_for_read(IIN,basedir,iproc,iregion,'cijkl_kernel')
  read(IIN) model_array
  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_write_model(model_array, basedir, iproc, iregion, model_name)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: model_array
  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  if (.not. is_dimension_set) then
    print *, 'sem_IO: undetermined dimension, involk sem_set_dimension() first!'
    stop
  end if

  call sem_open_datafile_for_write(IOUT,basedir,iproc,iregion,model_name)
  write(IOUT) model_array
  close(IOUT)

end subroutine

!//////////////////////////
subroutine sem_write_model_n(model_array, basedir, iproc, iregion, nmodel, model_names)

  integer, intent(in) :: nmodel
  real(kind=CUSTOM_REAL), dimension(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: model_array
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion

  integer :: i

  do i = 1, nmodel
    call sem_write_model(model_array(i,:,:,:,:), basedir, iproc, iregion, model_names(i))
  end do

end subroutine

!//////////////////////////
subroutine _zwgljd(z,w,np,alpha,beta)
!=======================================================================
!
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
!     -----------   weights associated with Jacobi polynomials of degree
!                   n = np-1.
!
!     Note : alpha and beta coefficients must be greater than -1.
!            Legendre polynomials are special case of Jacobi polynomials
!            just by setting alpha and beta to 0.
!
!=======================================================================

  implicit none

  double precision, parameter :: zero = 0.d0,one=1.d0,two=2.d0

  integer np
  double precision alpha,beta
  double precision z(np), w(np)

  integer n,nm1,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision alpg,betg
  double precision, external :: endw1,endw2

  p = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n   = np-1
  nm1 = n-1
  pd  = zero

  if (np <= 1) stop 'minimum number of Gauss-Lobatto points is 2'

! with spectral elements, use at least 3 points
  if (np <= 2) stop 'minimum number of Gauss-Lobatto points for the SEM is 3'

  if ((alpha <= -one) .or. (beta <= -one)) stop 'alpha and beta must be greater than -1'

  if (nm1 > 0) then
    alpg  = alpha+one
    betg  = beta+one
    call zwgjd(z(2),w(2),nm1,alpg,betg)
  endif

! start and end point at exactly -1 and 1
  z(1)  = - one
  z(np) =  one

! if number of points is odd, the middle abscissa is exactly zero
  if (mod(np,2) /= 0) z((np-1)/2+1) = zero

! weights
  do i=2,np-1
   w(i) = w(i)/(one-z(i)**2)
  enddo

  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
  w(1)  = endw1(n,alpha,beta)/(two*pd)
  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
  w(np) = endw2(n,alpha,beta)/(two*pd)

  end subroutine _zwgljd


end module sem_io
