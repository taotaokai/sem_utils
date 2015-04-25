module sem_io

!----------------------------------
! specification part
!----------------------------------

  use sem_constants

  implicit none

  private

  !---- public operations
  public :: sem_io_open_file_for_read
  public :: sem_io_open_file_for_write
  public :: sem_io_read_kernel_cijkl
  public :: sem_io_read_model
  public :: sem_io_read_modeln
  public :: sem_io_write_model
  public :: sem_io_write_modeln

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
subroutine sem_io_open_file_for_read(IIN,basedir,iproc,iregion,model_name)

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion, IIN

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(model_name)
  open(IIN,file=trim(filename),status='old', &
    form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'sem_io: Error open file for read: ', trim(filename)
    stop
  endif

end subroutine

!//////////////////////////
subroutine sem_io_open_file_for_write(IOUT,basedir,iproc,iregion,model_name)

  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion, IOUT

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  write(filename,"(a,'/proc',i6.6,'_reg',i1,'_',a,'.bin')") &
    trim(basedir), iproc, iregion, trim(model_name)
  open(IOUT,file=trim(filename),status='unknown', &
    form='unformatted',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'sem_io: Error open file for write: ', trim(filename)
    stop
  endif

end subroutine

!//////////////////////////
subroutine sem_io_read_model(model_array, basedir, iproc, iregion, model_name)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  integer :: ier

  call sem_io_open_file_for_read(IIN,basedir,iproc,iregion,model_name)
  read(IIN, iostat=ier) model_array
  if (ier /= 0) stop "sem_io_read_model:read error"

  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_io_read_modeln(model_array, basedir, iproc, iregion, nmodel, model_names)

  integer, intent(in) :: nmodel
  real(kind=CUSTOM_REAL), dimension(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion

  integer :: i

  do i = 1, nmodel
    call sem_io_read_model(model_array(i,:,:,:,:), basedir, iproc, iregion, model_names(i))
  end do

end subroutine

!//////////////////////////
subroutine sem_io_read_kernel_cijkl(model_array, basedir, iproc, iregion)

  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC), intent(out) :: model_array
  character(len=*), intent(in) :: basedir
  integer, intent(in) :: iproc, iregion

  call sem_io_open_datafile_for_read(IIN,basedir,iproc,iregion,'cijkl_kernel')
  read(IIN) model_array
  close(IIN)

end subroutine

!//////////////////////////
subroutine sem_io_write_model(model_array, basedir, iproc, iregion, model_name)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: model_array
  character(len=*), intent(in) :: basedir, model_name
  integer, intent(in) :: iproc, iregion

  call sem_io_open_datafile_for_write(IOUT,basedir,iproc,iregion,model_name)
  write(IOUT) model_array
  close(IOUT)

end subroutine

!//////////////////////////
subroutine sem_io_write_modeln(model_array, basedir, iproc, iregion, nmodel, model_names)

  integer, intent(in) :: nmodel
  real(kind=CUSTOM_REAL), dimension(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC), intent(in) :: model_array
  character(len=*), intent(in) :: basedir
  character(len=MAX_STRING_LEN), dimension(nmodel), intent(in) :: model_names
  integer, intent(in) :: iproc, iregion

  integer :: i

  do i = 1, nmodel
    call sem_io_write_model(model_array(i,:,:,:,:), basedir, iproc, iregion, model_names(i))
  end do

end subroutine


end module sem_io
