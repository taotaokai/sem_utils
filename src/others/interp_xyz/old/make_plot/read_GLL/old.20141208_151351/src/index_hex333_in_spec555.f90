! get index of anchor points(3x3x3 hex) in a 5x5x5 specfem

subroutine index_hex333_in_spec555(ixa,iya,iza)

  implicit none

! index of anchor points 
	integer, parameter :: NGNOD = 27
  integer, dimension(NGNOD) :: ixa,iya,iza

! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif

! corner nodes

  ixa(1) = 1
  iya(1) = 1
  iza(1) = 1

  ixa(2) = 5
  iya(2) = 1
  iza(2) = 1

  ixa(3) = 5
  iya(3) = 5
  iza(3) = 1

  ixa(4) = 1
  iya(4) = 5
  iza(4) = 1

  ixa(5) = 1
  iya(5) = 1
  iza(5) = 5

  ixa(6) = 5
  iya(6) = 1
  iza(6) = 5

  ixa(7) = 5
  iya(7) = 5
  iza(7) = 5

  ixa(8) = 1
  iya(8) = 5
  iza(8) = 5

! mid-edge nodes (nodes located in the middle of an edge)

  ixa(9) = 3
  iya(9) = 1
  iza(9) = 1

  ixa(10) = 5
  iya(10) = 3
  iza(10) = 1

  ixa(11) = 3
  iya(11) = 5
  iza(11) = 1

  ixa(12) = 1
  iya(12) = 3
  iza(12) = 1

  ixa(13) = 1
  iya(13) = 1
  iza(13) = 3

  ixa(14) = 5
  iya(14) = 1
  iza(14) = 3

  ixa(15) = 5
  iya(15) = 5
  iza(15) = 3

  ixa(16) = 1
  iya(16) = 5
  iza(16) = 3

  ixa(17) = 3
  iya(17) = 1
  iza(17) = 5

  ixa(18) = 5
  iya(18) = 3
  iza(18) = 5

  ixa(19) = 3
  iya(19) = 5
  iza(19) = 5

  ixa(20) = 1
  iya(20) = 3
  iza(20) = 5

! face center nodes (nodes located in the middle of a face)

  ixa(21) = 3
  iya(21) = 3
  iza(21) = 1

  ixa(22) = 3
  iya(22) = 1
  iza(22) = 3

  ixa(23) = 5
  iya(23) = 3
  iza(23) = 3

  ixa(24) = 3
  iya(24) = 5
  iza(24) = 3

  ixa(25) = 1
  iya(25) = 3
  iza(25) = 3

  ixa(26) = 3
  iya(26) = 3
  iza(26) = 5

! body center node (barycenter of the eight corners)

  ixa(27) = 3
  iya(27) = 3
  iza(27) = 3

end subroutine index_hex333_in_spec555
