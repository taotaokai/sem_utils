  ! shared variables
  integer :: i,j,k
  integer :: ipar
  integer :: iproc,ispec
  !integer :: igllx,iglly,igllz

  ! interpolation points
  integer :: ipt,npts
  real(RK), dimension(NDIM,NPTS_MAX) :: xyzi

  ! source mesh data
  integer :: nspec, nglob
  real(RK), dimension(NDIM,NGLOB_MAX) :: xyz
  real(RK), dimension(NDIM,NGNOD,NSPEC_MAX) :: xyza
  real(RK), dimension(NDIM,NSPEC_MAX) :: xyzc
  real(RK), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX,NPARS) :: model
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: ibool
  integer, dimension(NSPEC_MAX) :: idoubling

  ! grid location data for mesh 2
  real(RK), dimension(NDIM,NPTS_MAX) :: uvw
  integer, dimension(NPTS_MAX) :: eid
  real(RK), dimension(NPTS_MAX) :: res
  logical, dimension(NPTS_MAX) :: ins, fin

  ! locate_grid
  integer :: II
  integer :: n_sel
  real(RK), dimension(NSPEC_MAX) :: dist
  logical, dimension(NSPEC_MAX) :: idx_sel
  integer, dimension(NSPEC_MAX) :: ind_sel, ind_sort
  real(RK), dimension(NDIM) :: uvwi
  real(RK) :: resi
  logical :: inside, isdone

  ! interp_model
  real(RK) :: val
  real(RK), dimension(NGLLX) :: hxi, heta, hgamma
  real(RK), dimension(NPTS_MAX,NPARS) :: modeli

  ! write_model
  integer :: ier
  character(MAX_STRING_LEN) :: FMTOUT, basename, fn_model
