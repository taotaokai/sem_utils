module sem_io

!----------------------------------
! specification part
!----------------------------------

  use sem_constants, only: dp
  use sem_constants, only: CUSTOM_REAL
  use sem_constants, only: NGLLX, NGLLY, NGLLZ
  use sem_constants, only: MAX_STRING_LEN, IIN, IOUT

  implicit none

  private

  !---- public operations
  public :: sem_io_open_file_for_read
  public :: sem_io_open_file_for_write
  public :: sem_io_read_gll_file_1
  public :: sem_io_read_gll_file_n
  public :: sem_io_read_gll_file_cijkl
  public :: sem_io_write_gll_file_1
  public :: sem_io_write_gll_file_n

!----------------------------------
! implementation part
!----------------------------------
contains

!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_open_file_for_read(basedir, iproc, iregion, tag, fileid)
! open a GLL file, read only

  character(len=*), intent(in) :: basedir, tag 
  integer, intent(in) :: iproc, iregion, fileid

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(tag)

  open(unit=fileid, file=trim(filename), status='old', &
    form='unformatted',action='read',iostat=ier)

  if (ier /= 0) then
    write(*,'(a)') '[ERROR] sem_io_open_file_for_read: failed to open file ', &
                    trim(filename)
    stop
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_open_file_for_write(basedir, iproc, iregion, tag, fileid)
! open a GLL file to write

  character(len=*), intent(in) :: basedir, tag 
  integer, intent(in) :: iproc, iregion, fileid 

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename, "(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(tag)

  open(unit=fileid, file=trim(filename), status='unknown', &
    form='unformatted',action='write',iostat=ier)

  if (ier /= 0) then
    write(*,'(a)') '[ERROR] sem_io_open_file_for_write: failed to open file ', &
                    trim(filename)
    stop
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_read_gll_file_1(basedir, iproc, iregion, model_name, &
                                  model_gll)
!-read one GLL file into an array
!
!-inputs:
! (string) basedir,iproc,iregion,model_name: specify the GLL file
!
!-output:
! (real) model_gll(NGLLX,NGLLY,NGLLZ,NSPEC): GLL array

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  real(dp), intent(out) :: model_gll(:,:,:,:)

  integer :: ier, nspec
  real(kind=CUSTOM_REAL), allocatable :: dummy(:,:,:,:)

  nspec = size(model_gll,4)
  allocate(dummy(NGLLX,NGLLY,NGLLZ,nspec))

  call sem_io_open_file_for_read(basedir,iproc,iregion,model_name,IIN)

  read(IIN, iostat=ier) dummy

  if (ier /= 0) then
    stop "[ERROR] sem_io_read_gll_1: failed to read in gll file"
  endif

  close(IIN)

  model_gll = real(dummy, kind=dp)

end subroutine


!//////////////////////////
subroutine sem_io_read_gll_file_n(basedir, iproc, iregion, &
  model_names,nmodel, model_gll)
!-read n GLL files into one array
!
!-inputs:
! (string) basedir,iproc,iregion,model_names(nmodel): specify the GLL file
!
!-output:
! (real) model_gll(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC): GLL array

  integer, intent(in) :: nmodel
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion

  real(dp), intent(out) :: model_gll(:,:,:,:,:)

  integer :: imodel

  do imodel = 1, nmodel
    call sem_io_read_gll_file_1(basedir, iproc, iregion, model_names(imodel), &
                                model_gll(imodel,:,:,:,:))
  enddo

end subroutine


!//////////////////////////
subroutine sem_io_read_gll_file_cijkl(basedir, iproc, iregion, gll_cijkl)
!-read the special GLL file (c_ijkl) into one array
!
!-inputs:
! (string) basedir,iproc,iregion,model_names(nmodel): specify the GLL file
!
!-output:
! (real) model_gll(21,NGLLX,NGLLY,NGLLZ,NSPEC): GLL array

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  real(dp), intent(out) :: gll_cijkl(:,:,:,:,:)

  integer :: ier, nspec
  real(kind=CUSTOM_REAL), allocatable :: dummy(:,:,:,:,:)

  nspec = size(gll_cijkl, 5)
  allocate(dummy(21, NGLLX, NGLLY, NGLLZ, nspec))

  call sem_io_open_file_for_read(basedir, iproc, iregion, 'cijkl_kernel', IIN)

  read(IIN, iostat=ier) dummy 

  close(IIN)

  if (ier /= 0) then
    stop "[ERROR] sem_io_read_gll_file_ciikl: failed to read in gll file"
  endif

  gll_cijkl = real(dummy, kind=dp)

end subroutine


!//////////////////////////
subroutine sem_io_write_gll_file_1(basedir, iproc, iregion, &
  model_name, model_gll)
!-write out one GLL array onto disk
!
!-inputs:
! (string) basedir,iproc,iregion,model_name: specify the GLL file
!
!-output:
! (real) model_gll(:,:,:,:): GLL array

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion
  real(dp), intent(in) :: model_gll(:,:,:,:)

  integer :: ier

  call sem_io_open_file_for_write(basedir, iproc, iregion, model_name, IOUT)

  write(IOUT, iostat=ier) real(model_gll, kind=CUSTOM_REAL)

  close(IOUT)

  if (ier /= 0) then
    stop "[ERROR] sem_io_write_gll_file_1: failed to write out gll file"
  endif

end subroutine


!//////////////////////////
subroutine sem_io_write_gll_file_n(basedir, iproc, iregion, &
  model_names, nmodel, model_gll)
!-write out GLL array of multi-parameters into multiple files
!
!-inputs:
! (string) basedir,iproc,iregion,model_names(nmodel): specify the GLL files
!
!-output:
! (real) model_gll(nmodel,:,:,:,:): GLL array

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion
  integer, intent(in) :: nmodel
  character(len=MAX_STRING_LEN), intent(in) :: model_names(nmodel)

  real(dp), intent(in) :: model_gll(:,:,:,:,:)

  integer :: imodel

  do imodel = 1, nmodel
    call sem_io_write_gll_file_1(basedir, iproc, iregion, model_names(imodel), &
                                 model_gll(imodel,:,:,:,:))
  end do

end subroutine


end module sem_io
