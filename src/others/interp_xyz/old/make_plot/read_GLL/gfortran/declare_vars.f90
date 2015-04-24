  ! shared variables
  integer :: i1,j1,k1,i2
  integer :: iproc1,iproc2
  integer :: ipar
  integer :: ispec1,ispec2,iglob2
  integer :: igllx,iglly,igllz

  ! old mesh 1 data
  integer :: nspec, nglob
  real(RK), dimension(NDIM,NGLOB_MAX) :: xyz
  real(RK), dimension(NDIM,NGNOD,NSPEC_MAX) :: xyza
  real(RK), dimension(NDIM,NSPEC_MAX) :: xyzc
  real(RK), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX,NPARS) :: model
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: ibool

  ! location data in the mesh of each interpolation point
  real(RK), dimension(NDIM,NPOINT) :: uvw
  integer, dimension(NPOINT) :: eid
  real(RK), dimension(NPOINT) :: res
  logical, dimension(NPOINT) :: ins

  ! locate_point
  integer :: n_sel
  real(RK), dimension(NSPEC_MAX) :: dist
  logical, dimension(NSPEC_MAX) :: idx,idx_sel
  integer, dimension(NSPEC_MAX) :: ind_sel, ind_sort
  real(RK), dimension(NDIM) :: uvwi
  real(RK) :: resi
  logical :: inside

  ! interp_model
  real(RK) :: val
  real(RK), dimension(NGLLX) :: hxi, heta, hgamma

  ! write_model
  integer :: ier
  character(MAX_STRING_LEN) :: basename, fn_model
