! TIMESTAMP prints the current YMDHMS date as a time stamp.
! yyyy-mm-ddThh:mm:ss.sss

subroutine timestamp(out_date_string)

  use constants,only: MAX_STRING_LEN

  implicit none
  
  character(len=MAX_STRING_LEN), intent(out) :: out_date_string 

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

  write(out_date_string,'(i4,"-",i2,"-",i2,"T",i2,":",i2,":",i2,".",i3)') &
    y, m, d, h, n, s, ms

end subroutine