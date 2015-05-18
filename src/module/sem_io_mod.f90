module sem_io

!----------------------------------
! specification part
!----------------------------------

  use sem_constants, only: CUSTOM_REAL, MAX_STRING_LEN, IIN, IOUT

  implicit none

  private

  !---- public operations
  public :: sem_io_open_file_read
  public :: sem_io_open_file_write
  public :: sem_io_read_gll_1
  public :: sem_io_read_gll_n
  public :: sem_io_read_gll_cijkl
  public :: sem_io_write_gll_1
  public :: sem_io_write_gll_n

!----------------------------------
! implementation part
!----------------------------------
contains

!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_open_file_read(fileid, basedir, iproc, iregion, tag)
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
    write(*,'(a)') '[ERROR] sem_io_open_file_read: failed to open file ', &
                    trim(filename)
    stop
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_open_file_write(fileid, basedir, iproc, iregion, tag)
! open a GLL file to write

  character(len=*), intent(in) :: basedir, tag 
  integer, intent(in) :: iproc, iregion, fileid 

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename, "(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(tag)

  open(fileid, file=trim(filename), status='unknown', &
    form='unformatted',action='write',iostat=ier)

  if (ier /= 0) then
    write(*,'(a)') '[ERROR] sem_io_open_file_write: failed to open file ', &
                    trim(filename)
    stop
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_read_gll_1(basedir, iproc, iregion, model_name, &
                             dims, model_gll)
!-read one GLL file into an array
!
!-inputs:
! (string) basedir,iproc,iregion,model_name: specify the GLL file
! (integer) dims(4) = [NGLLX,NGLLY,NGLLZ,NSPEC]
!   specify the array dimension in the GLL file
!   (i.e. dimensions of GLL points of the SEM mesh)
!
!-output:
! (real) model_gll(NGLLX,NGLLY,NGLLZ,NSPEC): GLL array

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion, dims(4)

  real(kind=CUSTOM_REAL), intent(out) :: &
    model_gll(dims(1), dims(2), dims(3), dims(4))

  integer :: ier

  call sem_io_open_file_read(IIN,basedir,iproc,iregion,model_name)

  read(IIN, iostat=ier) model_gll

  if (ier /= 0) stop "[ERROR] sem_io_read_gll_1: failed to read in gll file"

  close(IIN)

end subroutine


!//////////////////////////
subroutine sem_io_read_gll_n(basedir, iproc, iregion, nmodel, model_names, &
                             dims, model_gll)
!-read n GLL files into one array
!
!-inputs:
! (string) basedir,iproc,iregion,model_names(nmodel): specify the GLL file
! (integer) dims(4) = [NGLLX,NGLLY,NGLLZ,NSPEC]
!   specify the array dimension in the GLL file
!   (i.e. dimensions of GLL points of the SEM mesh)
!
!-output:
! (real) model_gll(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC): GLL array

  integer, intent(in) :: nmodel
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion, dims(4)

  real(kind=CUSTOM_REAL), intent(out) :: &
    model_gll(nmodel, dims(1), dims(2), dims(3), dims(4))

  integer :: imodel

  do imodel = 1, nmodel
    call sem_io_read_gll_1(basedir, iproc, iregion, model_names(imodel) &
                           dims, model_gll(imodel,:,:,:,:))
  enddo

end subroutine


!//////////////////////////
subroutine sem_io_read_gll_cijkl(basedir, iproc, iregion, dims, gll_cijkl)
!-read the special GLL file (c_ijkl) into one array
!
!-inputs:
! (string) basedir,iproc,iregion,model_names(nmodel): specify the GLL file
! (integer) dims(4) = [NGLLX,NGLLY,NGLLZ,NSPEC]
!   specify the array dimension in the GLL file
!   (i.e. dimensions of GLL points of the SEM mesh)
!
!-output:
! (real) model_gll(21,NGLLX,NGLLY,NGLLZ,NSPEC): GLL array

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion, dims(4)

  real(kind=CUSTOM_REAL), intent(out) :: &
    gll_cijkl(21, dims(1), dims(2), dims(3), dims(4))

  integer :: ier

  call sem_io_open_file_read(IIN, basedir,iproc,iregion,'cijkl_kernel')

  read(IIN, iostat=ier) gll_cijkl

  close(IIN)

  if (ier /= 0) stop "[ERROR] sem_io_read_gll_ciikl: failed to read in gll file"

end subroutine


!//////////////////////////
subroutine sem_io_write_gll_1( &
  basedir, iproc, iregion, model_name, model_gll)
!-write out one GLL array onto disk
!
!-inputs:
! (string) basedir,iproc,iregion,model_name: specify the GLL file
!
!-output:
! (real) model_gll(:,:,:,:): GLL array

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion
  real(kind=CUSTOM_REAL), intent(in) :: model_gll(:,:,:,:)

  call sem_io_open_file_write(IOUT,basedir,iproc,iregion,model_name)

  write(IOUT, iostat=ier) model_gll

  close(IOUT)

  if (ier /= 0) stop "[ERROR] sem_io_write_gll_1: failed to write out gll file"

end subroutine


!//////////////////////////
subroutine sem_io_write_gll_n( &
  basedir, iproc, iregion, nmodel, model_names, model_gll)
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

  real(kind=CUSTOM_REAL), intent(in) :: model_gll(nmodel,:,:,:,:)

  integer :: imodel

  do imodel = 1, nmodel
    call sem_io_write_gll_1(basedir, iproc, iregion, model_names(imodel), &
                            model_gll(imodel,:,:,:,:))
  end do

end subroutine


end module sem_io
