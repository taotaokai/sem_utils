!------------------------------------------------------------
!    Implementation of a Heap Sort Routine
!
!    Input
!      n = Input
!         Length of arrays
!      X = Input/Output
!         Vector to be sorted
!         dimension(n)
!      Y = Output
!         Sorted Indices of vector X
!
!    Example:
!         X = [ 4 3 1 2 ] on Input
!         Y = [ 1 2 3 4 ] Computed Internally (in order)
!
!         X = [ 1 2 3 4 ] on Output
!         Y = [ 3 4 2 1 ] on Output
!------------------------------------------------------------

  subroutine heap_sort_local(N, X, Y)

  use sem_constants, only: CUSTOM_REAL

  implicit none

  integer, intent(in) :: N
  real(kind=CUSTOM_REAL), dimension(N), intent(inout) :: X
  integer, dimension(N), intent(out) :: Y

  ! local parameters
  real(kind=CUSTOM_REAL) :: tmp
  integer :: itmp
  integer :: i

  do i = 1,N
     Y(i) = i
  enddo

  ! checks if anything to do
  if (N < 2) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown(i, N)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    tmp = X(1)
    X(1) = X(i)
    X(i) = tmp
    itmp = Y(1)
    Y(1) = Y(i)
    Y(i) = itmp

    call heap_sort_siftdown(1, i - 1)
  enddo

!
!----
!

  contains

    subroutine heap_sort_siftdown(start, bottom)

    implicit none

    integer, intent(in) :: start, bottom

    ! local parameters
    integer :: i, j
    real(kind=CUSTOM_REAL) :: xtmp
    integer :: ytmp

    i = start
    xtmp = X(i)
    ytmp = Y(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
      if (j < bottom) then
        if (X(j) <= X(j+1)) j = j + 1
      endif

      ! checks if section already smaller than initial value
      if (X(j) < xtmp) exit

      X(i) = X(j)
      Y(i) = Y(j)
      i = j
      j = 2 * i
    enddo

    X(i) = xtmp
    Y(i) = ytmp

    end subroutine heap_sort_siftdown

  end subroutine heap_sort_local
