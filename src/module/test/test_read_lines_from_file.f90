program test

  integer, parameter :: MAX_STRING_LEN = 512

  character(len=MAX_STRING_LEN) :: filename

  integer :: nline 
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: slist

  integer :: i, ier
  character(len=MAX_STRING_LEN) :: dummy

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0) then
    filename = "text.lst"
  endif

  open(unit=IIN, file=trim(filename), status='old', iostat=ier)
  if (ier /= 0) then
     write(*,*) 'ERROR:text_utils:read_text_list: fail to open file: ', trim(filename)
     stop
  endif

  nline = 0 ! get number of lines
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    dummy = adjustl(dummy)
    if (len_trim(dummy) /= 0 .and. dummy(1:1) /= "#") then
      nline = nline + 1
      print *, "n=", nline
      print *, "line=", trim(dummy)
    endif
  enddo

! if (.not. allocated(slist)) then
!   allocate(slist(nline))
! elseif (size(slist) /= nline) then
!   deallocate(slist)
!   allocate(slist(nline))
! end if
! rewind(IIN)
! do i = 1, nline 
!   read(IIN,'(a)',iostat=ier) slist(i)
! enddo

  close(IIN)

end program