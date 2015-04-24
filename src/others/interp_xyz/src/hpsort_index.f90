!*****************************************************
!*  Index an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*      N     size of table RA                       *
!*      RA    table to be sorted                     *
!* OUTPUT:                                           *
!*      IX    index table sorted in ascending order  *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
SUBROUTINE HPSORT_INDEX(N,RA,IX)

  use constants, only: RK=>SIZE_REAL

  integer,intent(in) :: N
  real(RK),intent(in) :: RA(N)
  integer,intent(out) :: IX(N)

  L=N/2+1
  IR=N
  do I=1,N
    IX(I) = I
  end do

  if (N<2) return 

  do

    if(L>1)then
      L=L-1
      IIX=IX(L)
    else
      IIX=IX(IR)
      IX(IR)=IX(1)
      IR=IR-1
      if(IR.eq.1)then
        IX(1)=IIX
        return
      end if
    end if
  
    I=L
    J=L+L
  
    do while(J <= IR)
      if(J < IR)then
        if(RA(IX(J)) < RA(IX(J+1)))  J=J+1
      end if
      if(RA(IIX) < RA(IX(J)))then
        IX(I)=IX(J)
        I=J; J=J+J
      else
        J=IR+1
      end if
    end do
  
    IX(I)=IIX

  end do
END SUBROUTINE HPSORT_INDEX
