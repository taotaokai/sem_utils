module utils

!----------------------------------
! specification part
!----------------------------------

  use constants,only: MAX_STRING_LEN,IIN,IOUT

  implicit none

private

  !---- public

  ! subroutines
  public :: read_text_list
  public :: timestamp
  public :: delimit_string 

!----------------------------------
! implementation part
!----------------------------------
contains

!//////////////////////////
! read text file line by line
subroutine read_text_list(filename, slist, nline)

  character(len=MAX_STRING_LEN), intent(in) :: filename

  integer, intent(out) :: nline 
  character(len=MAX_STRING_LEN), dimension(:), allocatable, intent(out) :: slist

  integer :: i, ier
  character(len=MAX_STRING_LEN) :: dummy

  open(unit=IIN, file=trim(filename), status='old', iostat=ier)
  if (ier /= 0) then
     write(*,*) 'Error open file: ', trim(filename)
     stop
  endif
  nline = 0 ! get number of lines
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    nline = nline + 1
  enddo
  if (.not. allocated(slist)) then
    allocate(slist(nline))
  elseif (size(slist) /= nline) then
    deallocate(slist)
    allocate(slist(nline))
  end if
  rewind(IIN)
  do i = 1, nline 
    read(IIN,'(a)',iostat=ier) slist(i)
  enddo
  close(IIN)

end subroutine

!//////////////////////////
! time stamp: yyyy-mm-ddThh:mm:ss.sss
subroutine timestamp(out_str)

  implicit none

  character(len=MAX_STRING_LEN), intent(out) :: out_str

  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  integer :: values(8)
  integer y, m, d, h, n, s, ms

  call date_and_time(date, time, zone, values)
  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  ms = values(8)

  write(out_str,'(i4,"-",i2,"-",i2,"T",i2,":",i2,":",i2,".",i3)') &
    y, m, d, h, n, s, ms
end subroutine !timestamp

!//////////////////////////
subroutine delimit_string(source_str, delimiter, str_array, n)

  implicit none

  character(len=*), intent(in) :: source_str
  character(len=*), intent(in) :: delimiter
  character(len=MAX_STRING_LEN), intent(inout), allocatable :: str_array(:)
  integer, intent(out) :: n

  ! local parameters
  integer :: i
  character(len=MAX_STRING_LEN) :: dummy_str 

  ! get number of fields in the string
  n = 0
  do
    call _strtok_(source_str, delimiter, dummy)
    if ((dummy(1) /= char(0) .and. (len_trim(dummy)/=0)) then
      n = n + 1
    endif
  enddo

  ! allocate array
  if (allocated(str_array)) deallocate(str_array)
  allocate(str_array(n))

  ! split string into array
  do i = 1, n
    call _strtok_(source_str, delimiter, dummy)
    if ((dummy(1) /= char(0) .and. (len_trim(dummy)/=0)) then
      str_arrayn(i) = dummy
    endif
  enddo

end subroutine delimit_string

!//////////////////////////
! The following utility function was modified from http://fortranwiki.org/fortran/show/strtok
subroutine _strtok_(source_string, delimiter, token)
! @(#) Tokenize a string in a similar manner to C routine strtok(3c).
!
! Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!         and the delimiter list used to tokenize SOURCE_STRING in delimiter.
!
!         then, if the returned value is not equal to char(0), keep calling until it is
!         with SOURCE_STRING set to char(0).
!
!        STRTOK will return a token on each call until the entire line is processed,
!        which it signals by returning char(0).
!
! Input:  source_string =   Source string to tokenize.
!         delimiter    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
! Output: strtok()
!
! LIMITATIONS:
! can not be called with a different string until current string is totally processed, even from different procedures

  !     PARAMETERS:
  character(len=*), intent(in)  :: source_string
  character(len=*), intent(in)  :: delimiter
  character(len=MAX_STRING_LEN), intent(out) :: token

  !     SAVED VALUES:
  character(len=MAX_STRING_LEN),save :: saved_string
  integer,save :: isaved_start  ! points to beginning of unprocessed data
  integer,save :: isource_len   ! length of original input string

  !     LOCAL VALUES:
  integer :: ibegin        ! beginning of token to return
  integer :: ifinish       ! end of token to return

  ! initialize stored copy of input string and pointer into input string on first call
  if (source_string(1:1) /= char(0)) then
    isaved_start = 1                 ! beginning of unprocessed data
    saved_string = source_string     ! save input string from first call in series
    isource_len = LEN(saved_string)  ! length of input string from first call
  endif

  token = ''
  ibegin = isaved_start

  ! sets first index ibegin to beginning of (next) token
  do while (.true.)
    if ( (ibegin <= isource_len) .and. (index(delimiter,saved_string(ibegin:ibegin)) /= 0)) then
      ! delimiter is encountered, starts with next index (next token)
      ibegin = ibegin + 1
    else
      ! exits do-loop
      exit
    endif
  enddo

  if (ibegin > isource_len) then
    token = char(0)
    return
  endif

  ! sets second index ifinish to end of token (including delimiter)
  ifinish = ibegin

  do while (.true.)
    if ((ifinish <= isource_len) .and.  (index(delimiter,saved_string(ifinish:ifinish)) == 0)) then
      ! delimiter is not encountered yet, increases finish index
      ifinish = ifinish + 1
    else
      ! exits do-loop
      exit
    endif
  enddo

  ! sets token string
  !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
  token = saved_string(ibegin:ifinish-1)
  isaved_start = ifinish

end subroutine _strtok_

end module
