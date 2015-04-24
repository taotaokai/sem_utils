! test formated output

program test

!integer, parameter :: IOUT = 0

integer :: a = 1
real :: x = 2.345

!open(unit=IOUT, file='out.txt', status='unknown', iostat=ier)

do i = 1,3
  write(*, '($,I1,2x,E15.7)') a,x
enddo

write(*, '()')

do i = 1,3
  write(*, '($,I1,2x,E15.7)') a,x
enddo

!close(IOUT)

end program
