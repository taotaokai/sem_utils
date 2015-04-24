module sem_IO

!----------------------------------
! specification part
!----------------------------------

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ, &
    IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE

  implicit none

private

  ! array dimension
  integer :: NSPEC, NGLOB
  logical :: flag_set_dimension=.false.
  ! initial flag
  logical :: flag_init_mesh = .false.

  ! model structure
  type, public :: model_iso
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
      mu, lamda_2mu, rho
  end type

  ! mesh structure
  type, public :: mesh
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyz
    integer, dimension(:,:,:,:), allocatable :: ibool
    integer, dimension(:), allocatable :: idoubling
    logical, dimension(:), allocatable :: ispec_is_tiso 
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  end type mesh

  !---- public
  ! dimensions 
  public :: NGLLX,NGLLY,NGLLZ,NSPEC,NGLOB

  ! subroutines
  public :: sem_set_dimension
  public :: sem_open_datafile_for_read
  public :: sem_open_datafile_for_write
  public :: sem_init_mesh
  public :: sem_read_mesh
  public :: sem_read_kernel_cijkl
  public :: sem_read_model
  public :: sem_write_model
  public :: sem_read_model_iso
  public :: sem_write_model_iso

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

  flag_set_dimension = .true.

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

  type (mesh), intent(inout) :: mesh_data

  if (.not. flag_set_dimension) then
    print *, 'sem_IO: undetermined dimension, invoke sem_set_dimension() first!'
    stop
  end if

  if (allocated(mesh_data%xyz)) then
    return
  end if

  allocate(mesh_data%xyz(3,NGLOB))
  allocate(mesh_data%ibool(NGLLX,NGLLY,NGLLZ,NSPEC), &
           mesh_data%idoubling(NSPEC), &
           mesh_data%ispec_is_tiso(NSPEC))
  allocate(mesh_data%xix(NGLLX,NGLLY,NGLLZ,NSPEC), &
           mesh_data%xiy(NGLLX,NGLLY,NGLLZ,NSPEC), &
           mesh_data%xiz(NGLLX,NGLLY,NGLLZ,NSPEC), &
          mesh_data%etax(NGLLX,NGLLY,NGLLZ,NSPEC), &
          mesh_data%etay(NGLLX,NGLLY,NGLLZ,NSPEC), &
          mesh_data%etaz(NGLLX,NGLLY,NGLLZ,NSPEC), &
        mesh_data%gammax(NGLLX,NGLLY,NGLLZ,NSPEC), &
        mesh_data%gammay(NGLLX,NGLLY,NGLLZ,NSPEC), &
        mesh_data%gammaz(NGLLX,NGLLY,NGLLZ,NSPEC))

  flag_init_mesh = .true.

end subroutine

!//////////////////////////
subroutine sem_read_mesh(mesh_data, basedir, iproc, iregion)

  type (mesh), intent(inout) :: mesh_data
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  integer :: ival

  if (.not. flag_set_dimension) then
    print *, 'sem_IO: undetermined dimension, involk sem_set_dimension() first!'
    stop
  end if
  if (.not. flag_init_mesh) then
    print *, 'sem_IO: type(mesh) not initialized, involk sem_init_mesh() first!'
    stop
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

  read(IIN) mesh_data%xyz(1,:)
  read(IIN) mesh_data%xyz(2,:)
  read(IIN) mesh_data%xyz(3,:)

  read(IIN) mesh_data%ibool(NGLLX,NGLLY,NGLLZ,NSPEC)
  read(IIN) mesh_data%idoubling(1:NSPEC)
  read(IIN) mesh_data%ispec_is_tiso(1:NSPEC)

  read(IIN) mesh_data%xix
  read(IIN) mesh_data%xiy
  read(IIN) mesh_data%xiz
  read(IIN) mesh_data%etax
  read(IIN) mesh_data%etay
  read(IIN) mesh_data%etaz
  read(IIN) mesh_data%gammax
  read(IIN) mesh_data%gammay
  read(IIN) mesh_data%gammaz

  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_read_model(model_array, basedir, iproc, iregion, model_name)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  integer :: ier

  if (.not. flag_set_dimension) then
    print *, 'sem_IO: undetermined dimension, involk sem_set_dimension() first!'
    stop
  end if

  call sem_open_datafile_for_read(IIN,basedir,iproc,iregion,model_name)
  read(IIN, iostat=ier) model_array
  if (ier /= 0) stop "sem_read_model:read error"

  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_read_model_iso(model_data, basedir, iproc, iregion)

  type(model_iso), intent(inout) :: model_data
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  if (.not. allocated(model_data%mu)) then
    allocate(model_data%mu(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(model_data%lamda_2mu(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(model_data%rho(NGLLX,NGLLY,NGLLZ,NSPEC))
  end if

  call sem_read_model(model_data%lamda_2mu, basedir, iproc, iregion, "lamda_2mu")
  call sem_read_model(model_data%mu, basedir, iproc, iregion, "mu")
  call sem_read_model(model_data%rho, basedir, iproc, iregion, "rho")

end subroutine

!//////////////////////////
subroutine sem_read_kernel_cijkl(model_array, basedir, iproc, iregion)

  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  if (.not. flag_set_dimension) then
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

  if (.not. flag_set_dimension) then
    print *, 'sem_IO: undetermined dimension, involk sem_set_dimension() first!'
    stop
  end if

  call sem_open_datafile_for_write(IOUT,basedir,iproc,iregion,model_name)
  write(IOUT) model_array
  close(IOUT)

end subroutine

!//////////////////////////
subroutine sem_write_model_iso(model_data, basedir, iproc, iregion)

  type(model_iso), intent(inout) :: model_data
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  call sem_write_model(model_data%mu, basedir, iproc, iregion, "mu")
  call sem_write_model(model_data%lamda_2mu, basedir, iproc, iregion, "lamda_2mu")
  call sem_write_model(model_data%rho, basedir, iproc, iregion, "rho")

end subroutine

end module sem_IO
