program test_string

  implicit none

  integer, parameter :: MAX_STRING_LEN = 512

  character(len=MAX_STRING_LEN) :: source_string, dummy, arg
  character(len=*), parameter :: delimiter = ","

  integer :: n

  call get_command_argument(1, source_string)
  if (len_trim(source_string) == 0) then
    source_string = "this, is, a, #test"
  endif

  n = 0
  dummy = adjustl(strtok(source_string, delimiter))
  do
    if (dummy == char(0)) then
      exit
    else
      if (len_trim(dummy) /= 0 .and. dummy(1:1) /= "#") then
        n = n + 1
        print *, "n=", n
        print *, "dummy=", trim(dummy), "|END"
      endif
    endif
    dummy = adjustl(strtok(char(0), delimiter))
  enddo

contains

CHARACTER(len=MAX_STRING_LEN) FUNCTION strtok (source_string, delimiters)
!     @(#) Tokenize a string in a similar manner to C routine strtok(3c). 
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in DELIMITERS.
!
!             then, if the returned value is not equal to CHAR(0), keep calling until it is
!             with SOURCE_STRING set to CHAR(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning CHAR(0). 
!
!     Input:  source_string =   Source string to tokenize. 
!             delimiters    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!     Output: strtok()
!
!     LIMITATIONS:
!     can not be called with a different string until current string is totally processed, even from different procedures
!     input string length limited to set size
!     function returns fixed 255 character length
!     length of returned string not given

!     PARAMETERS:
      CHARACTER(len=*),intent(in)  :: source_string
      CHARACTER(len=*),intent(in)  :: delimiters

!     SAVED VALUES:
      CHARACTER(len=255),save :: saved_string
      INTEGER,save :: isaved_start  ! points to beginning of unprocessed data
      INTEGER,save :: isource_len   ! length of original input string

!     LOCAL VALUES:
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
         IF ( (ibegin .LE. isource_len) .AND. (INDEX(delimiters,saved_string(ibegin:ibegin)) .NE. 0)) THEN
             ibegin = ibegin + 1
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
         IF ((ifinish .LE. isource_len) .AND.  (INDEX(delimiters,saved_string(ifinish:ifinish)) .EQ. 0)) THEN
             ifinish = ifinish + 1
         ELSE
             EXIT
         ENDIF
      ENDDO

      !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
      strtok = saved_string(ibegin:ifinish-1)
      isaved_start = ifinish

END FUNCTION strtok

end program