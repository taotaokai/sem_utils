subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle - create slice of SEM model along great circle"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_tags> "
  print '(a)', "    <lat0> <lon0> <azimuth> "
  print '(a)', "    <theta0> <theta1> <ntheta> "
  print '(a)', "    <radius0> <radius1> <nradius> "
  print '(a)', "    <flag_ellipticity> <out_file> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  create a verical slice of SEM model along specified great circle" 
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc: number of mesh slices"
  print '(a)', "  (string) mesh_dir: directory containing proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir: directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0,lon0: origin point on the great circle (degrees)"
  print '(a)', "  (float) azimuth: shooting azimuth of the great circle (degrees)"
  print '(a)', "  (float) theta0, theta1: begin/end angles along great circle (degrees)"
  print '(a)', "  (integer) ntheta: number of grids along the great cicle"
  print '(a)', "  (float) radius0,radius1: begin/end radius (km)"
  print '(a)', "  (integer) nradius: number of grids along the radius"
  print '(a)', "  (int) flag_ellipticity: (0/1) 0: spherical; 1: WGS84 ellipsoid"
  print '(a)', "  (string) out_file: output file name (GMT netcdf format)"
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  1. Example of plot: bash scripts using GMT5"
  print '(a)', "    ~~~~~~~~~~~~~~~~"
  print '(a)', "    #!/bin/bash "
  print '(a)', "    ps=<model_tag>.ps "
  print '(a)', "    R=<theta0>/<theta1>/<radius0>/<radius1>"
  print '(a)', "    gmt grdreformat <out_file>?<model_tag> <model_tag>.grd "
  print '(a)', "    gmt grd2cpt <model_tag>.grd -Cseis -R$R -Z -D > cpt "
  print '(a)', "    gmt grdimage <model_tag>.grd -R$R -JPa6iz -Ccpt \"
  print '(a)', "      -BNSEW -Bxa10f5 -Bya100f50 > $ps "
  print '(a)', "    gs $ps" 
  print '(a)', "    ~~~~~~~~~~~~~~~~"
  print '(a)', ""
  print '(a)', "  2. To run in parallel: "
  print '(a)', "    mpirun -n <nproc_mpi> <program> <arguments>"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_vertical_slice
  
  use sem_constants
  use sem_parallel
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic
  use netcdf

  implicit none

  !===== declare variables

  !-- command line args
  integer, parameter :: nargs = 15
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_file
  real(dp) :: lat0, lon0, azimuth, theta0, theta1, radius0, radius1
  integer :: nproc, ntheta, nradius
  integer :: flag_ellipticity

  !-- local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, ier, iproc

  !-- model names
  integer :: imodel, nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- interpolation points 
  integer :: ipoint, npoint
  real(dp) :: v0(3), v1(3), vnorth(3), veast(3), v_axis(3), vr(3), rotmat(3,3)
  real(dp) :: dtheta, dradius
  real(dp), allocatable :: theta(:), radius(:), lat(:), lon(:)
  integer :: itheta, iradius, idx
  real(dp), allocatable :: xyz(:,:)
  integer, allocatable :: idoubling(:)

  !-- sem location 
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec, iglob
  real(dp) :: dist
  type(sem_mesh_location), allocatable :: location_1slice(:)
  integer, parameter :: nnearest = 10
  real(dp) :: typical_size, max_search_dist, max_misloc
  real(dp), allocatable :: misloc_myrank(:), misloc_gather(:), misloc_copy(:)
  integer, allocatable :: stat_myrank(:), stat_gather(:), stat_copy(:)
  !-- model interpolation
  real(sp), parameter :: FILLVALUE_sp = huge(0.0_sp)
  real(dp), allocatable :: model_gll(:,:,:,:,:)
  real(dp), allocatable :: model_interp_myrank(:,:)
  real(dp), allocatable :: model_interp_gather(:,:)
  real(dp), allocatable :: model_interp_copy(:,:)

  !-- mpi 
  integer :: myrank, nrank, iworker
  ! mpi_send/recv
  integer, parameter :: MPI_TAG_stat = 10
  integer, parameter :: MPI_TAG_misloc = 11
  integer, parameter :: MPI_TAG_model_interp = 12

  !-- netcdf data structure
  integer :: ncid
  ! dimensions: theta, radius
  integer, parameter :: NDIMS = 2
  integer :: theta_dimid, radius_dimid, dimids(NDIMS)
  ! coordinates: theta(:), radius(:), lat(:), lon(:)
  integer :: theta_varid, radius_varid, lat_varid, lon_varid
  character(len=*), parameter :: UNITS = "units"
  character(len=*), parameter :: theta_name = "theta"
  character(len=*), parameter :: theta_units = "degree"
  character(len=*), parameter :: radius_name = "radius"
  character(len=*), parameter :: radius_units = "km"
  character(len=*), parameter :: lat_name = "latitude"
  character(len=*), parameter :: lat_units = "degree"
  character(len=*), parameter :: lon_name = "longitude"
  character(len=*), parameter :: lon_units = "degree"
  ! slice geometry
  integer :: lat0_varid, lon0_varid, az0_varid
  character(len=*), parameter :: lat0_name = "latitude_orig"
  character(len=*), parameter :: lat0_units = "degree"
  character(len=*), parameter :: lon0_name = "longitude_orig"
  character(len=*), parameter :: lon0_units = "degree"
  character(len=*), parameter :: az0_name = "azimuth_orig"
  character(len=*), parameter :: az0_units = "degree"
  ! data: x,y,z(ntheta, nradius)
  integer :: x_varid, y_varid, z_varid
  real(sp), dimension(:,:), allocatable :: x_var, y_var, z_var
  character(len=*), parameter :: x_name = "x"
  character(len=*), parameter :: x_units = "R_EARTH"
  character(len=*), parameter :: y_name = "y"
  character(len=*), parameter :: y_units = "R_EARTH"
  character(len=*), parameter :: z_name = "z"
  character(len=*), parameter :: z_units = "R_EARTH"
  ! data: model values (units is not determined)
  integer, allocatable :: model_varids(:)
  real(sp), allocatable :: model_var(:,:)
  ! data: location status, mislocation
  integer :: stat_varid, misloc_varid
  integer, allocatable :: stat_var(:,:)
  real(sp), allocatable :: misloc_var(:,:)
  character(len=*), parameter :: stat_name = "location_status"
  character(len=*), parameter :: misloc_name = "misloc"

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_slice_gcircle: check your input arguments."
      call abort_mpi()
      stop
    endif
  endif
   
  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), *) lat0
  read(args(6), *) lon0
  read(args(7), *) azimuth
  read(args(8), *) theta0
  read(args(9), *) theta1 
  read(args(10), *) ntheta
  read(args(11), *) radius0
  read(args(12), *) radius1
  read(args(13), *) nradius
  read(args(14), *) flag_ellipticity
  read(args(15), '(a)') out_file

  call synchronize_all()

  !===== generate mesh grid of the cross-section

  ! convert units, non-dimensionalize as in SPECFEM3D_GLOBE 
  lat0 = lat0 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  azimuth = azimuth * DEGREES_TO_RADIANS

  theta0 = theta0 * DEGREES_TO_RADIANS
  theta1 = theta1 * DEGREES_TO_RADIANS

  radius0 = radius0 / R_EARTH_KM
  radius1 = radius1 / R_EARTH_KM

  ! unit radial vector v0 from earth's center to the origin point on the great circle
  if (flag_ellipticity /= 0) then
    call geographic_lla2ecef(lat0, lon0, 0.0_dp, v0(1), v0(2), v0(3))
    v0 = v0 / sqrt(sum(v0**2))
  else ! spherical earth
    v0(1) = cos(lat0)*cos(lon0)
    v0(2) = cos(lat0)*sin(lon0)
    v0(3) = sin(lat0)
  endif

  ! unit direction vector v1 along the shooting azimuth of the great circle
  vnorth = [ - sin(lat0) * cos(lon0), &
             - sin(lat0) * sin(lon0), &
               cos(lat0) ]
  veast = [ - sin(lon0), cos(lon0), 0.d0 ]

  v1 = cos(azimuth) * vnorth + sin(azimuth) * veast

  ! rotation axis = v0 cross-product v1
  v_axis(1) = v0(2)*v1(3) - v0(3)*v1(2)
  v_axis(2) = v0(3)*v1(1) - v0(1)*v1(3)
  v_axis(3) = v0(1)*v1(2) - v0(2)*v1(1)
  v_axis = v_axis / sqrt(sum(v_axis**2))

  ! get grid intervals
  dtheta = (theta1 - theta0) / (ntheta - 1)
  dradius = (radius1 - radius0) / (nradius - 1)

  ! mesh grid xyz(1:3, radius*theta)
  allocate(theta(ntheta), radius(nradius), lat(ntheta), lon(ntheta))
  theta = [ (theta0 + i * dtheta, i=0,ntheta-1) ]
  radius = [ (radius0 + i * dradius, i=0,nradius-1 ) ]

  npoint = nradius * ntheta
  allocate(xyz(3, npoint))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    vr = matmul(rotmat, v0)

    if (flag_ellipticity == 0) then
      lat(itheta) = atan2(vr(3), (vr(1)**2+vr(2)**2)**0.5)
      lon(itheta) = atan2(vr(2), vr(1))
    else
      call geographic_ecef2ll_zeroalt(vr(1), vr(2), vr(3), lat(itheta), lon(itheta))
    endif

    idx = nradius * (itheta - 1)
    do iradius = 1, nradius
      xyz(:, idx+iradius) = radius(iradius) * vr
    enddo
  enddo

  ! layer ID 
  allocate(idoubling(npoint))
  idoubling = IFLAG_DUMMY

  call synchronize_all()

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, delimiter, model_names, nmodel)

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)
  endif

  call synchronize_all()

  !===== locate xyz in each mesh slice

  !-- initialize variables
  allocate(location_1slice(npoint))

  allocate(stat_myrank(npoint), misloc_myrank(npoint))
  stat_myrank = -1
  misloc_myrank = HUGE(1.0_dp)

  allocate(model_interp_myrank(nmodel, npoint))
  model_interp_myrank = FILLVALUE_sp

  ! typical element size at surface
  typical_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                     ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                 * DEGREES_TO_RADIANS * R_UNIT_SPHERE

  max_search_dist = 5.0 * typical_size

  max_misloc = typical_size / 4.0

  !-- loop mesh slices for myrank
  loop_proc: do iproc = myrank, (nproc - 1), nrank

    print *, "# iproc=", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! test if the interpolation points and mesh are separated apart
    ! only use element centers
    dist = huge(0.0_dp)
    do ispec = 1, nspec
      iglob = mesh_data%ibool(MIDX, MIDY, MIDZ, ispec)
      dist = min(dist, sqrt( minval( &
        (xyz(1,:) - mesh_data%xyz_glob(1,iglob))**2 + &
        (xyz(2,:) - mesh_data%xyz_glob(2,iglob))**2 + &
        (xyz(3,:) - mesh_data%xyz_glob(3,iglob))**2)))
    enddo
    if (dist > max_search_dist) then
      cycle
    endif

    ! read model
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(nmodel, NGLLX, NGLLY, NGLLZ, nspec))

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
      model_names, nmodel, model_gll)

    ! locate points in this mesh slice
    call sem_mesh_locate_kdtree2(mesh_data, npoint, xyz, idoubling, &
      nnearest, max_search_dist, max_misloc, location_1slice)

    print *, "# iproc=", iproc, &
      " number of located points: ", count(location_1slice%stat>=0)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! safety check
      if (stat_myrank(ipoint) == 1 .and. &
          location_1slice(ipoint)%stat == 1 ) then
        print *, "[WARN] iproc=", iproc, " ipoint=", ipoint, &
          " This point is located inside more than one element!", &
          " Only use the first located element."
        cycle
      endif

      ! for point located inside one element but not happened before 
      ! or closer to one element than located before
      if (location_1slice(ipoint)%stat == 1 &
          .or. &
          (location_1slice(ipoint)%stat == 0 .and. &
           location_1slice(ipoint)%misloc < misloc_myrank(ipoint)) ) &
      then

        ! interpolate model
        do imodel = 1, nmodel
          model_interp_myrank(imodel, ipoint) = &
            sum( location_1slice(ipoint)%lagrange * &
                 model_gll(imodel, :, :, :, location_1slice(ipoint)%eid) )
        enddo

        stat_myrank(ipoint)   = location_1slice(ipoint)%stat
        misloc_myrank(ipoint) = location_1slice(ipoint)%misloc

      endif

    enddo ! ipoint

  enddo loop_proc

  !-- gather all location results from work processors
  call synchronize_all()

  if (myrank /= 0) then

    ! send stat/misloc/model_interp on myrank to MASTER 

    call send_i(stat_myrank, npoint, 0, MPI_TAG_stat)
    call send_dp(misloc_myrank, npoint, 0, MPI_TAG_misloc)
    call send2_dp(model_interp_myrank, nmodel*npoint, 0, MPI_TAG_model_interp)

  else ! gather data on MASTER processor (0)

    allocate(stat_gather(npoint), misloc_gather(npoint), &
      model_interp_gather(nmodel, npoint))

    allocate(stat_copy(npoint), misloc_copy(npoint), &
      model_interp_copy(nmodel, npoint))

    stat_gather = stat_myrank 
    misloc_gather = misloc_myrank 
    model_interp_gather = model_interp_myrank 

    do iworker = 1, (nrank -1)

      ! receive stat/misloc/model_interp from each worker processes 
      call recv_i(stat_copy, npoint, iworker , MPI_TAG_stat)
      call recv_dp(misloc_copy, npoint, iworker, MPI_TAG_misloc)
      call recv2_dp(model_interp_copy, nmodel*npoint, iworker, &
        MPI_TAG_model_interp)

      ! gather data
      do ipoint = 1, npoint

        ! safety check
        if (stat_gather(ipoint) == 1) then
          if (stat_copy(ipoint) == 1) then
            print *, "[WARN] ipoint=", ipoint
            print *, "------ this point is located inside more than one element!"
            print *, "------ some problem may occur."
            print *, "------ only use the first located element."
          endif
          cycle
        endif

        ! for point located inside one element but not happened before 
        ! or closer to one element than located before
        if (stat_copy(ipoint) == 1 &
            .or. &
            (stat_copy(ipoint) == 0 .and. &
             misloc_copy(ipoint) < misloc_gather(ipoint)) ) &
        then
          stat_gather(ipoint) = stat_copy(ipoint)
          misloc_gather(ipoint) = misloc_copy(ipoint)
          model_interp_gather(:, ipoint) = model_interp_copy(:,ipoint) 
        endif

      enddo ! ipoint

    enddo ! iworker

    ! convert misloc relative to typical element size
    where (stat_gather /= -1)
      misloc_gather = misloc_gather / typical_size
    endwhere

  endif ! myrank /= 0

  !===== output netcdf file on the MASTER processor

  if (myrank == 0) then

    !-- create the file
    call check( nf90_create(out_file, NF90_CLOBBER, ncid))

    !-- define dimensions.
    call check( nf90_def_dim(ncid, theta_name, ntheta, theta_dimid) )
    call check( nf90_def_dim(ncid, radius_name, nradius, radius_dimid) )

    !-- define coordinate variables.
    call check( nf90_def_var(ncid, theta_name, NF90_DOUBLE, theta_dimid, theta_varid) )
    call check( nf90_def_var(ncid, radius_name, NF90_DOUBLE, radius_dimid, radius_varid) )
    call check( nf90_def_var(ncid, lat_name, NF90_DOUBLE, theta_dimid, lat_varid) )
    call check( nf90_def_var(ncid, lon_name, NF90_DOUBLE, theta_dimid, lon_varid) )
    ! coordinate units
    call check( nf90_put_att(ncid, theta_varid, UNITS, theta_units) )
    call check( nf90_put_att(ncid, radius_varid, UNITS, radius_units) )
    call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
    call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )

    !-- define slice geometry variables
    call check( nf90_def_var(ncid, lat0_name, NF90_DOUBLE, lat0_varid) )
    call check( nf90_put_att(ncid, lat0_varid, UNITS, lat0_units) )
    call check( nf90_def_var(ncid, lon0_name, NF90_DOUBLE, lon0_varid) )
    call check( nf90_put_att(ncid, lon0_varid, UNITS, lon0_units) )
    call check( nf90_def_var(ncid, az0_name, NF90_DOUBLE, az0_varid) )
    call check( nf90_put_att(ncid, az0_varid, UNITS, az0_units) )

    !-- define data variables.
    dimids = [theta_dimid, radius_dimid]
    ! x,y,z(ntheta,nradius)
    call check( nf90_def_var(ncid, x_name, NF90_REAL, dimids, x_varid) )
    call check( nf90_put_att(ncid, x_varid, UNITS, x_units) )
    call check( nf90_def_var(ncid, y_name, NF90_REAL, dimids, y_varid) )
    call check( nf90_put_att(ncid, y_varid, UNITS, y_units) )
    call check( nf90_def_var(ncid, z_name, NF90_REAL, dimids, z_varid) )
    call check( nf90_put_att(ncid, z_varid, UNITS, z_units) )
    ! model values(ntheta,nradius)
    allocate(model_varids(nmodel))
    do imodel = 1, nmodel
      call check( nf90_def_var(ncid, trim(model_names(imodel)), NF90_REAL & 
                               , dimids, model_varids(imodel)) )
      call check( nf90_put_att(ncid, model_varids(imodel) &
                               , "_FillValue", FILLVALUE_sp) )
    enddo
    ! location status and mislocation
    call check( nf90_def_var(ncid, stat_name, NF90_INT, dimids, stat_varid) )
    call check( nf90_def_var(ncid, misloc_name, NF90_REAL, dimids, misloc_varid) )

    ! end define mode.
    call check( nf90_enddef(ncid) )

    ! write coordinate variables
    call check( nf90_put_var(ncid, theta_varid, theta * RADIANS_TO_DEGREES) )
    call check( nf90_put_var(ncid, radius_varid, radius * R_EARTH_KM) )
    call check( nf90_put_var(ncid, lat_varid, lat * RADIANS_TO_DEGREES) )
    call check( nf90_put_var(ncid, lon_varid, lon * RADIANS_TO_DEGREES) )

    ! write slice geometry variables
    call check( nf90_put_var(ncid, lat0_varid, lat0 * RADIANS_TO_DEGREES) )
    call check( nf90_put_var(ncid, lon0_varid, lon0 * RADIANS_TO_DEGREES) )
    call check( nf90_put_var(ncid, az0_varid, azimuth * RADIANS_TO_DEGREES) )

    ! write data: x,y,z
    allocate(x_var(ntheta, nradius), y_var(ntheta, nradius), z_var(ntheta, nradius))
    do itheta = 1, ntheta
      idx = nradius * (itheta - 1)
      do iradius = 1, nradius
        x_var(itheta, iradius) = real(xyz(1, idx + iradius), kind=sp)
        y_var(itheta, iradius) = real(xyz(2, idx + iradius), kind=sp)
        z_var(itheta, iradius) = real(xyz(3, idx + iradius), kind=sp)
      enddo
    enddo
    call check( nf90_put_var(ncid, x_varid, x_var) )
    call check( nf90_put_var(ncid, y_varid, y_var) )
    call check( nf90_put_var(ncid, z_varid, z_var) )

    ! write data: interpolated model values
    allocate(model_var(ntheta, nradius))
    do imodel = 1, nmodel

      ! reshape model variable into a NDIMS matrix
      do itheta = 1, ntheta
        idx = nradius * (itheta - 1)
        do iradius = 1, nradius
          model_var(itheta, iradius) = &
            real(model_interp_gather(imodel, idx + iradius), kind=sp)
        enddo
      enddo

      call check( nf90_put_var(ncid, model_varids(imodel), model_var) )
    enddo

    ! write data: status of location and mislocation
    ! reshape variable into an NDIMS matrix
    allocate(stat_var(ntheta, nradius), misloc_var(ntheta, nradius))
    do itheta = 1, ntheta
      idx = nradius * (itheta - 1)
      do iradius = 1, nradius
        stat_var(itheta, iradius) = stat_gather(idx + iradius)
        misloc_var(itheta, iradius) = real(misloc_gather(idx + iradius), kind=sp)
      enddo
    enddo

    call check( nf90_put_var(ncid, stat_varid, stat_var) )
    call check( nf90_put_var(ncid, misloc_varid, misloc_var) )

    ! close file 
    call check( nf90_close(ncid) )

  endif ! myrank == 0

  !===== Finalize MPI

  call synchronize_all()
  call finalize_mpi()


!********************
contains

subroutine check(status)
  integer, intent (in) :: status
  
  if(status /= nf90_noerr) then 
    print "('[ERROR] ',a)", trim(nf90_strerror(status))
    stop
  endif
end subroutine check 


end program
