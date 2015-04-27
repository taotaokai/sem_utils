module sem_io

!----------------------------------
! specification part
!----------------------------------

  use sem_constants

  implicit none

  private

  !---- public operations
  public :: sem_io_open_proc_file_for_read
  public :: sem_io_open_proc_file_for_write
  public :: sem_io_read_kernel_cijkl
  public :: sem_io_read_gll_model
  public :: sem_io_read_gll_modeln
  public :: sem_io_write_gll_model
  public :: sem_io_write_gll_modeln

!----------------------------------
! implementation part
!----------------------------------
contains

!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_open_proc_file_for_read(fileid, basedir, iproc, iregion, tag)

  character(len=*), intent(in) :: basedir, tag 
  integer, intent(in) :: iproc, iregion, fileid 

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(tag)

  open(unit=fileid, file=trim(filename), status='old', &
    form='unformatted',action='read',iostat=ier)

  if (ier /= 0) then
    write(*,*) 'ERROR:sem_io_open_proc_file_for_read: failed to open file: ', trim(filename)
    stop
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_open_proc_file_for_write(fileid, basedir, iproc, iregion, tag)

  character(len=*), intent(in) :: basedir, tag 
  integer, intent(in) :: iproc, iregion, fileid 

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename, "(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(tag)

  open(fileid, file=trim(filename), status='unknown', &
    form='unformatted',action='write',iostat=ier)

  if (ier /= 0) then
    write(*,*) 'ERROR:sem_io_open_proc_file_for_write: failed to open file: ', trim(filename)
    stop
  endif

end subroutine


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_io_read_gll_model(basedir, iproc, iregion, model_name, model_gll)

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  real(kind=CUSTOM_REAL), intent(out) :: model_gll(NGLLX,NGLLY,NGLLZ,NSPEC)

  integer :: ier

  call sem_io_open_proc_file_for_read(IIN,basedir,iproc,iregion,model_name)
  read(IIN, iostat=ier) model_gll
  if (ier /= 0) stop "ERROR: sem_io_read_gll_model: failed to read gll model"

  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_io_read_gll_modeln(basedir, iproc, iregion, nmodel, model_names &
                                  , model_gll)

  integer, intent(in) :: nmodel
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion

  real(kind=CUSTOM_REAL), intent(out) :: model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel)

  integer :: imodel

  do imodel = 1, nmodel
    call sem_io_read_gll_model(basedir, iproc, iregion, model_names(imodel) &
                               , model_gll(:,:,:,:,imodel))
  end do

end subroutine


!//////////////////////////
subroutine sem_io_read_kernel_cijkl(basedir, iproc, iregion, kernel_cijkl)

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  real(kind=CUSTOM_REAL), intent(out) :: kernel_cijkl(21,NGLLX,NGLLY,NGLLZ,NSPEC)

  call sem_io_open_proc_file_for_read(IIN,basedir,iproc,iregion,'cijkl_kernel')
  read(IIN) kernel_cijkl
  close(IIN)

end subroutine


!//////////////////////////
subroutine sem_io_write_gll_model(basedir, iproc, iregion, model_name, model_gll)

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion
  real(kind=CUSTOM_REAL), intent(in) :: model_gll(NGLLX,NGLLY,NGLLZ,NSPEC)

  call sem_io_open_proc_file_for_write(IOUT,basedir,iproc,iregion,model_name)
  write(IOUT) model_gll
  close(IOUT)

end subroutine


!//////////////////////////
subroutine sem_io_write_gll_modeln( &
  basedir, iproc, iregion, nmodel, model_names, model_gll)

  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion
  integer, intent(in) :: nmodel
  character(len=MAX_STRING_LEN), intent(in) :: model_names(nmodel)
  real(kind=CUSTOM_REAL), intent(in) :: model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel)

  integer :: imodel

  do imodel = 1, nmodel
    call sem_io_write_gll_model(basedir, iproc, iregion, model_names(imodel) &
                                , model_gll(:,:,:,:,imodel))
  end do

end subroutine


end module sem_io
