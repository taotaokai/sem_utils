! The following utility function was made freely available by the Fortran Wiki:
! http://fortranwiki.org/fortran/show/strtok
!
character(len=255) function strtok (source_string, delimiters)

!     @(#) Tokenize a string in a similar manner to C routine strtok(3c).
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in DELIMITERS.
!
!             then, if the returned value is not equal to char(0), keep calling until it is
!             with SOURCE_STRING set to char(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning char(0).
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
      character(len=*),intent(in)  :: source_string
      character(len=*),intent(in)  :: delimiters

!     SAVED VALUES:
      character(len=255),save :: saved_string
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

      ibegin = isaved_start

      do
         if ( (ibegin <= isource_len) .AND. (index(delimiters,saved_string(ibegin:ibegin)) /= 0)) then
             ibegin = ibegin + 1
         else
             exit
         endif
      enddo

      if (ibegin > isource_len) then
          strtok = char(0)
          RETURN
      endif

      ifinish = ibegin

      do
         if ((ifinish <= isource_len) .AND.  (index(delimiters,saved_string(ifinish:ifinish)) == 0)) then
             ifinish = ifinish + 1
         else
             exit
         endif
      enddo

      !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
      strtok = saved_string(ibegin:ifinish-1)
      isaved_start = ifinish

end function strtok

