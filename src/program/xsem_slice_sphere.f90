subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slice_sphere - create spherical slice of SEM model at a constant depth"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slice_shpere \ " 
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_tags> "
  print '(a)', "    <lat0> <lat1> <nlat> "
  print '(a)', "    <lon0> <lon1> <nlon> "
  print '(a)', "    <depth> <out_file> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  create a spherical slice of SEM model at a constant depth" 
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc*_reg1_solver_data.bin"
  print '(a)', "  (int) nproc: number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0, lat1:  begin/end latitudes (degrees) "
  print '(a)', "    valid range: [-90, 90] "
  print '(a)', "  (integer) nlat:  number of latitudinal grids"
  print '(a)', "  (float) lon0, lon1:  begin/end longitudes (degrees)"
  print '(a)', "    valid range: [-180, 180] "
  print '(a)', "  (integer) nlon:  number of longitudinal grids"
  print '(a)', "  (float) depth:  depth(0: r=6371) of the spherical slice (km)"
  print '(a)', "  (string) out_file:  output file name (netcdf format)"
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  1. Example of plot: bash scripts using GMT5"
  print '(a)', "    ~~~~~~~~~~~~~~~~"
  print '(a)', "    #!/bin/bash "
  print '(a)', "    ps=<model_tag>.ps "
  print '(a)', "    R=<lon0>/<lon1>/<lat0>/<lat1>"
  print '(a)', "    gmt grdreformat <out_file>?<model_tag> <model_tag>.grd "
  print '(a)', "    gmt grd2cpt <model_tag>.grd -Cseis -R$R -Z -D > cpt "
  print '(a)', "    gmt grdimage <model_tag>.grd -R$R -JM6i -Ccpt \"
  print '(a)', "      -BnSeW -Bxa10f5 -Bya10f5 > $ps "
  print '(a)', "    gs $ps"
  print '(a)', "    ~~~~~~~~~~~~~~~~"
  print '(a)', ""
  print '(a)', "  2. To run in parallel: "
  print '(a)', "    mpirun -n <nproc_mpi> <program> <arguments>"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_slice_sphere
  
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
  integer, parameter :: nargs = 12
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_file
  real(dp) :: lat0, lat1, lon0, lon1, depth 
  integer :: nlat, nlon, nproc

  !-- local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, ier, iproc

  !-- model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- interpolation points 
  integer :: ipoint, npoint
  real(dp) :: dlon, dlat, radius
  real(dp) :: vr(3)
  real(dp), allocatable :: lon(:), lat(:)
  integer :: ilon, ilat, idx
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
  ! dimensions: lon, lat
  integer, parameter :: NDIMS = 2
  integer :: lon_dimid, lat_dimid, dimids(NDIMS)
  ! coordinates: lon(:), lat(:)
  integer :: lon_varid, lat_varid
  character(len=*), parameter :: UNITS = "units"
  character(len=*), parameter :: lon_name = "longitude"
  character(len=*), parameter :: lon_units = "degrees_east"
  character(len=*), parameter :: lat_name = "latitude"
  character(len=*), parameter :: lat_units = "degrees_north"
  ! slice geometry
  integer :: depth_varid
  character(len=*), parameter :: depth_name = "depth"
  character(len=*), parameter :: depth_units = "km"
  ! define data variables
  ! this is defined from model_tags
  ! data units is not determined
  integer, allocatable :: model_varids(:)
  real(sp), allocatable :: model_var(:,:)
  ! define data: location status, mislocation
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
      print *, "[ERROR] xsem_slice_sphere: check your input arguments."
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
  read(args(6), *) lat1
  read(args(7), *) nlat
  read(args(8), *) lon0
  read(args(9), *) lon1
  read(args(10), *) nlon
  read(args(11), *) depth 
  read(args(12), '(a)') out_file

  call synchronize_all()

  !-- check validity of parameters

  !-90 <= lat0 < lat1 <=90
  if ( .not.(lat0 >= -90.d0 .and. lat0 <= 90.d0 .and. &
             lat1 >= -90.d0 .and. lat1 <= 90.d0 .and. &
             lat0 < lat1) ) then
    print *, "[ERROR] wrong inputs: lat0=", lat0, " lat1=", lat1
    call abort_mpi()
    stop
  endif

  !-180 <= lon0,lon1 <=180
  if ( .not.(lon0 >= -180.d0 .and. lat0 <= 180.d0 .and. &
             lon1 >= -180.d0 .and. lat1 <= 180.d0) ) then
    print *, "[ERROR] wrong ranges of input lat/lon"
    call abort_mpi()
    stop
  endif
  ! make sure the great circle path follows the Earth's rotation direction
  ! from lon0 to lon1
  if (lon0 > lon1) then
    print '(a)', "[WARN] lon0 > lon1, so 360 is added to lon1 "
    lon1 = lon1 + 360.d0
  endif

  !===== generate mesh grid of the cross-section

  ! convert units, non-dimensionalize as in SPECFEM3D_GLOBE 
  lat0 = lat0 * DEGREES_TO_RADIANS
  lat1 = lat1 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  lon1 = lon1 * DEGREES_TO_RADIANS
  radius = (6371.0 - depth) / R_EARTH_KM

  ! get grid intervals
  dlon = (lon1 - lon0) / (nlon - 1)
  dlat = (lat1 - lat0) / (nlat - 1)

  ! mesh grid xyz(1:3, radius*theta)
  allocate(lon(nlon), lat(nlat))
  lon = [ (lon0 + i * dlon, i = 0, nlon-1) ]
  lat = [ (lat0 + i * dlat, i = 0, nlat-1) ]

  npoint = nlon * nlat
  allocate(xyz(3,npoint))

  do ilat = 1, nlat
    do ilon = 1, nlon
      idx = ilon + nlon * (ilat - 1)
      call geographic_lla2ecef(lat(ilat), lon(ilon), 0.0d0, &
                                    vr(1), vr(2), vr(3)) 
      vr = vr / sqrt(sum(vr**2))
      xyz(:,idx) = vr * radius
    enddo
  enddo

  ! layer ID 
  allocate(idoubling(npoint))
  idoubling = IFLAG_DUMMY

  call synchronize_all()

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

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
    ! only use mesh element centers
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

    ! create the file
    call check( nf90_create(out_file, NF90_CLOBBER, ncid))

    ! define the dimensions.
    call check( nf90_def_dim(ncid, lon_name, nlon, lon_dimid) )
    call check( nf90_def_dim(ncid, lat_name, nlat, lat_dimid) )

    ! define the coordinate variables.
    call check( nf90_def_var(ncid, lon_name, NF90_DOUBLE, lon_dimid, lon_varid) )
    call check( nf90_def_var(ncid, lat_name, NF90_DOUBLE, lat_dimid, lat_varid) )
    ! coordinate units
    call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )
    call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )

    !-- define slice geometry variable
    call check( nf90_def_var(ncid, depth_name, NF90_DOUBLE, depth_varid) )
    call check( nf90_put_att(ncid, depth_varid, UNITS, depth_units) )

    ! define data variables: model values
    dimids = [lon_dimid, lat_dimid]
    allocate(model_varids(nmodel))
    do imodel = 1, nmodel
      call check( nf90_def_var(ncid, trim(model_names(imodel)), NF90_REAL & 
                               , dimids, model_varids(imodel)) )
      call check( nf90_put_att(ncid, model_varids(imodel) &
                               , "_FillValue", FILLVALUE_sp) )
    enddo
    ! define data: location status and mislocation
    call check( nf90_def_var(ncid, stat_name, NF90_INT, dimids, stat_varid) )
    call check( nf90_def_var(ncid, misloc_name, NF90_FLOAT, dimids, misloc_varid) )

    ! end define mode.
    call check( nf90_enddef(ncid) )

    ! write the coordinate variable data.
    lon = lon * RADIANS_TO_DEGREES
    lat = lat * RADIANS_TO_DEGREES
    ! change longitude back to range [-180,180]
    where (lon > 360.d0) lon = lon - 360.d0

    call check( nf90_put_var(ncid, lon_varid, lon) )
    call check( nf90_put_var(ncid, lat_varid, lat) )

    ! write slice geometry variable
    call check( nf90_put_var(ncid, depth_varid, depth) )

    ! write data: interpolated model values
    allocate(model_var(nlon, nlat))
    do imodel = 1, nmodel

      ! reshape model variable into a NDIMS matrix
      do ilat = 1, nlat
        idx = nlon * (ilat - 1)
        do ilon = 1, nlon
          model_var(ilon, ilat) = &
            real(model_interp_gather(imodel, idx + ilon), kind=sp)
        enddo
      enddo

      call check( nf90_put_var(ncid, model_varids(imodel), model_var) )
    enddo

    ! write data: status of location and mislocation
    ! reshape variable into an NDIMS matrix
    allocate(stat_var(nlon, nlat), misloc_var(nlon, nlat))
    do ilat = 1, nlat
      idx = nlon*(ilat - 1)
      do ilon = 1, nlon
        stat_var(ilon, ilat) = stat_gather(idx + ilon)
        misloc_var(ilon, ilat) = real(misloc_gather(idx + ilon), kind=sp)
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
