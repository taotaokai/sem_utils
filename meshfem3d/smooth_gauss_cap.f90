subroutine smooth_gauss_cap(ngll_target, xyz_gll_target, &
                            ngll_contrib, xyz_gll_contrib, vol_gll_contrib, &
                            nmodel, model_gll_contrib, sigmaH, sigmaV, &
                            weight_model_out, weight_out)

implicit none

! input/output
integer, intent(in) :: ngll_target, ngll_contrib, nmodel

! integer intent(hide),depend(xyz_gll_target) :: ngll_target = shape(xyz_gll_target)(2)
! integer intent(hide),depend(xyz_gll_contrib) :: ngll_contrib = shape(xyz_gll_contrib)(2)
! integer intent(hide),depend(model_gll_contrib) :: nmodel = shape(model_gll_contrib)(1)

real(kind=8), intent(in) :: xyz_gll_target(3,ngll_target)
real(kind=8), intent(in) :: xyz_gll_contrib(3,ngll_contrib)
real(kind=8), intent(in) :: model_gll_contrib(nmodel,ngll_contrib)
real(kind=8), intent(in) :: vol_gll_contrib(ngll_contrib)
real(kind=8), intent(in) :: sigmaH(ngll_target), sigmaV(ngll_target)

real(kind=8), intent(out) :: weight_model_out(nmodel,ngll_target)
real(kind=8), intent(out) :: weight_out(ngll_target)

! local variables
integer :: igll_target, igll_contrib
real(kind=8) :: r_target, r_gll_contrib(ngll_contrib)
real(kind=8) :: iprod, sigma2_V, dist2_theta, dist2_radial
real(kind=8) :: sigma2_theta, weight

! initialization
weight_model_out = 0.0
weight_out = 0.0

r_gll_contrib = sum(xyz_gll_contrib**2, 1)**0.5

do igll_target = 1, ngll_target

  r_target = sum(xyz_gll_target(:,igll_target)**2)**0.5

  sigma2_theta = (sigmaH(igll_target)/r_target)**2
  sigma2_V = sigmaV(igll_target)**2

  ! replace negative/zero sigma with a very small positive value
  !if (sigma2_theta <= 0.d0)  sigma2_theta = tiny(sigma2_theta)
  !if (sigma2_V <= 0.d0)  sigma2_V = tiny(sigma2_V)

  do igll_contrib = 1, ngll_contrib

    dist2_radial = (r_gll_contrib(igll_contrib) - r_target)**2

    iprod = sum(xyz_gll_contrib(:,igll_contrib)*xyz_gll_target(:,igll_target)) &
            / r_gll_contrib(igll_contrib) / r_target

    if (iprod > 1) iprod = 1.0
    dist2_theta = acos(iprod)**2

    iprod = -0.5*(dist2_radial/sigma2_V + dist2_theta/sigma2_theta)
    ! skip too small exponent in case of underflow(??)
    if (iprod > -10) then
      weight = exp(iprod) * vol_gll_contrib(igll_contrib)
      weight_model_out(:,igll_target) = weight_model_out(:,igll_target) + weight*model_gll_contrib(:,igll_contrib)
      weight_out(igll_target) = weight_out(igll_target) + weight
    endif

  enddo

enddo

end subroutine