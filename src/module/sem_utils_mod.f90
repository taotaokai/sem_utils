module sem_utils

!----------------------------------
! specification part
!----------------------------------

  use sem_constants,only: MAX_STRING_LEN, IIN, IOUT

  implicit none

private

  !---- public subroutines
  public :: sem_utils_read_line
  public :: sem_utils_delimit_string
  public :: sem_utils_timestamp

!----------------------------------
! implementation part
!----------------------------------
contains

!///////////////////////////////////////////////////////////////////////////////
subroutine sem_utils_read_line(filename, lines, n)
! read lines from a text file and output into an array of strings

  character(len=MAX_STRING_LEN), intent(in) :: filename

  integer, intent(out) :: n
  character(len=MAX_STRING_LEN), allocatable, intent(out) :: lines(:)

  integer :: ier
  character(len=MAX_STRING_LEN) :: dummy

  open(unit=IIN, file=trim(filename), status='old', iostat=ier)
  if (ier /= 0) then
     write(*,*) 'ERROR:file_readline: fail to open file: ', trim(filename)
     stop
  endif

  n = 0 ! get number of lines
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    dummy = adjustl(dummy)
    if (len_trim(dummy) /= 0 .and. dummy(1:1) /= "#") then
      n = n + 1
    endif
  enddo

  if (allocated(lines)) deallocate(lines)
  allocate(lines(n))

  rewind(IIN)

  n = 0 ! get number of lines
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    dummy = adjustl(dummy)
    if (len_trim(dummy) /= 0 .and. dummy(1:1) /= "#") then
      n = n + 1
      lines(n) = trim(dummy)
    endif
  enddo

end subroutine


!///////////////////////////////////////////////////////////////////////////////
function sem_utils_timestamp() result(string_out)
! timestamp: yyyy-mm-ddThh:mm:ss.sss

  character(len=MAX_STRING_LEN) :: string_out

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

  write(string_out,'(i0.4,"-",i0.2,"-",i0.2,"T",i0.2,":",i0.2,":",i0.2,".",i0.3)') &
        y, m, d, h, n, s, ms

end function 


!///////////////////////////////////////////////////////////////////////////////
subroutine sem_utils_delimit_string(source_string, delimiter, split_array, n)
! split delimited string into an array of strings

  character(len=*), intent(in) :: source_string
  character(len=*), intent(in) :: delimiter
  character(len=MAX_STRING_LEN), intent(inout), allocatable :: split_array(:)
  integer, intent(out) :: n

  ! local parameters
  character(len=MAX_STRING_LEN) :: dummy

  ! get number of fields in the string
  n = 0
  dummy = adjustl(strtok(source_string, delimiter))
  do
    if (dummy == char(0)) then
      exit
    else
      if (len_trim(dummy) /= 0) then
        n = n + 1
      endif
    endif
    dummy = adjustl(strtok(char(0), delimiter))
  enddo

  ! allocate array
  if (allocated(split_array)) deallocate(split_array)
  allocate(split_array(n))

  ! split string into array
  n = 0
  dummy = adjustl(strtok(source_string, delimiter))
  do
    if (dummy == char(0)) then
      exit
    else
      if (len_trim(dummy) /= 0) then
        n = n + 1
        split_array(n) = trim(dummy)
      endif
    endif
    dummy = adjustl(strtok(char(0), delimiter))
  enddo

end subroutine


!///////////////////////////////////////////////////////////////////////////////
character(len=MAX_STRING_LEN) function strtok (source_string, delimiters)
!-@(#) Tokenize a string in a similar manner to C routine strtok(3c). 
!
!-Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!         and the delimiter list used to tokenize SOURCE_STRING in DELIMITERS.
!
!         then, if the returned value is not equal to CHAR(0), keep calling it
!         with SOURCE_STRING set to CHAR(0).
!
!        STRTOK will return a token on each call until the entire line is processed,
!        which it signals by returning CHAR(0). 
!
!-Input:  source_string =   Source string to tokenize. 
!         delimiters    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!-Output: strtok()
!
!-LIMITATIONS:
! can not be called with a different string until current string is totally processed, even from different procedures
! input string length limited to set size
! function returns fixed 255 character length
! length of returned string not given

  ! PARAMETERS:
  CHARACTER(len=*),intent(in)  :: source_string
  CHARACTER(len=*),intent(in)  :: delimiters

  ! SAVED VALUES:
  CHARACTER(len=MAX_STRING_LEN),save :: saved_string
  INTEGER,save :: isaved_start  ! points to beginning of unprocessed data
  INTEGER,save :: isource_len   ! length of original input string

  ! LOCAL VALUES:
  INTEGER :: ibegin        ! beginning of token to return
  INTEGER :: ifinish       ! end of token to return

  ! initialize stored copy of input string and pointer into input string on first call
  IF (source_string(1:1) .NE. CHAR(0)) THEN
    isaved_start = 1                 ! beginning of unprocessed data
    saved_string = source_string     ! save input string from first call in series
    isource_len = LEN(saved_string)  ! length of input string from first call
  ENDIF

  ibegin = isaved_start

  DO
    IF (ibegin .LE. isource_len) THEN 
      IF (INDEX(delimiters,saved_string(ibegin:ibegin)) .NE. 0) THEN
        ibegin = ibegin + 1
      ELSE
        EXIT
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO

  IF (ibegin .GT. isource_len) THEN
    strtok = CHAR(0)
    RETURN
  ENDIF

  ifinish = ibegin

  DO
    IF (ifinish .LE. isource_len) THEN
      IF (INDEX(delimiters,saved_string(ifinish:ifinish)) .EQ. 0) THEN
        ifinish = ifinish + 1
      ELSE
        EXIT
      ENDIF
    ELSE
        EXIT
    ENDIF
  ENDDO

  !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
  strtok = saved_string(ibegin:ifinish-1)
  isaved_start = ifinish

end function strtok


end module
