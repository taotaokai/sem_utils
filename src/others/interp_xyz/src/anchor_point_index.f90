! get the index of anchor points in the GLL points of spectral element
! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif

! currently anchor points are defined on 3x3x3 hex

subroutine anchor_point_index(iax,iay,iaz)

  use constants,only: NGLLX,NGLLY,NGLLZ,NGNOD

  implicit none

  !==========================================================
  ! specification of subroutine arguments 
  !==========================================================

  ! index of anchor points
  integer :: ia
  integer, dimension(NGNOD), intent(out) :: iax,iay,iaz

  !==========================================================
  ! declare local variables 
  !==========================================================

  ! topology of the control points of the surface element
  integer, dimension(NGNOD) :: iaddx, iaddy, iaddz

  !==========================================================
  ! set up variables 
  !==========================================================

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddz)

  !==========================================================
  ! real work 
  !==========================================================

  do ia = 1,NGNOD

    if (iaddx(ia) == 0) then
      iax(ia) = 1
    else if (iaddx(ia) == 1) then
      iax(ia) = (NGLLX+1)/2
    else if (iaddx(ia) == 2) then
      iax(ia) = NGLLX
    else
      stop 'incorrect value of iaddx'
    endif

    if (iaddy(ia) == 0) then
      iay(ia) = 1
    else if (iaddy(ia) == 1) then
      iay(ia) = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay(ia) = NGLLY
    else
      stop 'incorrect value of iaddy'
    endif

    if (iaddz(ia) == 0) then
      iaz(ia) = 1
    else if (iaddz(ia) == 1) then
      iaz(ia) = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
      iaz(ia) = NGLLZ
    else
      stop 'incorrect value of iaddz'
    endif

  end do ! do ia = 1,NGNOD

end subroutine anchor_point_index
