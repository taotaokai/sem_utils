module test 

  type :: mesh(NGLOB)
    integer, len :: NGLOB
    integer, dimension(NGLOB) :: x, y, z
  end type

!  type mesh2
!    integer :: m, n, l
!  end type

end