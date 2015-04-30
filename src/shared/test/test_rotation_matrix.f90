program test

  use geographic

  implicit none

  real(8) :: lat0, lon0, lat1, lon1, r0, r1
  integer :: i, ntheta, nradius, npoint, idx, itheta, iradius

  real(8) :: v0(3), v1(3), dtheta, dradius, rotmat(3,3), v_axis(3)
  real(8), allocatable :: theta(:), radius(:), xyz(:,:)

  lat0 = 10.0
  lon0 = 10.0
  lat1 = 30.0
  lon1 = 90.0
  ntheta = 100

  r0 = (6371.0 - 1000.0) * 1000.0
  r1 = 6371.0 * 1000.0
  nradius = 100

  !integer :: i, j
  !real(8) :: v_axis(3), theta, R(3,3), v0(3), v1(3)

  !v_axis = [0.0, 0.0, 1.0]

  !theta = PI / 2.0

  !call rotation_matrix(v_axis, theta, R)

  !v0 = [1.0, 0.0, 1.0]

  !do i = 1,3
  !  do j = 1,3
  !    print *, i, j, R(i,j)
  !  enddo
  !enddo

  !v1 = matmul(R, v0)

  !print *, "v1=", v1

  !---- generate mesh grid of the cross-section
  call geodetic2ecef(lat0,lon0,0.0d0,v0)
  call geodetic2ecef(lat1,lon1,0.0d0,v1)

  print *, v0
  print *, v1

  ! normalize v0,v1 to uint vector 
  v0 = v0 / sqrt(sum(v0**2))
  v1 = v1 / sqrt(sum(v1**2))

  ! rotation axis = v0 cross-product v1
  v_axis(1) = v0(2)*v1(3) - v0(3)*v1(2)
  v_axis(2) = v0(3)*v1(1) - v0(1)*v1(3)
  v_axis(3) = v0(1)*v1(2) - v0(2)*v1(1)

  ! get gird intervals
  dtheta = acos(sum(v0*v1)) / (ntheta - 1)
  dradius = (r1 - r0) / (nradius - 1)

  ! get mesh grid
  allocate(theta(ntheta), radius(nradius))
  theta = (/(i * dtheta, i=0,ntheta-1 )/)
  radius = (/(r0 + i*dradius, i=0,nradius-1 )/)

  npoint = nradius * ntheta
  allocate(xyz(3,npoint))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    do iradius = 1, nradius
      idx = iradius + nradius*(itheta - 1)
      xyz(:,idx) = radius(iradius) * matmul(rotmat, v0)
      print *, theta(itheta), radius(iradius), xyz(:,idx)
    enddo
  enddo

end program