subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slice_sphere - create spherical slice of SEM model at a constant radius"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slice_shpere <mesh_dir> <model_dir> <model_tags> "
  print '(a)', "                    <lat0> <lat1> <nlat> "
  print '(a)', "                    <lon0> <lon1> <nlon> "
  print '(a)', "                    <radius> <out_file> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  create a spherical slice of SEM model at a constant radius" 
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0, lat1:  begin/end latitudes (degrees) "
  print '(a)', "  (integer) nlat:  number of latitudinal grids"
  print '(a)', "  (float) lon0, lon1:  begin/end longitudes (degrees)"
  print '(a)', "  (integer) nlon:  number of longitudinal grids"
  print '(a)', "  (float) raidus:  radius of the spherical slice (km)"
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

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_slice_sphere
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic
  use netcdf

  implicit none

  !===== declare variables

  integer, parameter :: dp = 8

  !-- command line args
  integer, parameter :: nargs = 11
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_file
  real(dp) :: lat0, lat1, lon0, lon1, radius
  integer :: nlat, nlon

  !-- local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  real(CUSTOM_REAL), parameter :: FILLVALUE = huge(0._CUSTOM_REAL)
  integer :: i, iproc, ier

  !-- model names
  integer :: imodel, nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- interpolation points 
  integer :: ipoint, npoint
  real(dp) :: dlon, dlat
  real(dp) :: vr(3)
  real(dp), allocatable :: lon(:), lat(:)
  integer :: ilon, ilat, idx
  real(CUSTOM_REAL), allocatable :: xyz(:,:)

  !-- sem interpolation
  type (mesh) :: mesh_data
  real(CUSTOM_REAL), allocatable :: model_gll(:,:,:,:,:), model_interp(:,:)
  real(CUSTOM_REAL), allocatable :: uvw(:,:), hlagrange(:,:,:,:)
  real(CUSTOM_REAL), allocatable :: misloc(:), misloc_final(:)
  integer, allocatable :: elem_ind(:), statloc(:)
  integer, allocatable :: statloc_final(:), elem_ind_final(:)
  real(CUSTOM_REAL) :: typical_size

  !-- netcdf data structure
  integer :: ncid
  ! define dimensions: theta, radius
  integer, parameter :: NDIMS = 2
  integer :: lon_dimid, lat_dimid, dimids(NDIMS)
  ! define coordinate variables/units
  integer :: lon_varid, lat_varid
  character(len=*), parameter :: UNITS = "units"
  character(len=*), parameter :: lon_name = "longitude"
  character(len=*), parameter :: lon_units = "degrees_east"
  character(len=*), parameter :: lat_name = "latitude"
  character(len=*), parameter :: lat_units = "degrees_north"
  ! define data variables
  ! this is defined from model_tags
  ! data units is not determined
  integer, allocatable :: model_varids(:)
  real(CUSTOM_REAL), allocatable :: model_var(:,:)
  ! define data: location status, mislocation
  integer :: statloc_varid, misloc_varid
  integer, allocatable :: statloc_var(:,:)
  real(CUSTOM_REAL), allocatable :: misloc_var(:,:)
  character(len=*), parameter :: statloc_name = "statloc"
  character(len=*), parameter :: misloc_name = "misloc"

  !===== read command line arguments

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (len_trim(args(i)) == 0) then
      call selfdoc()
      stop "ERROR: check your input arguments!"
    endif
  enddo

  read(args(1), '(a)') mesh_dir
  read(args(2), '(a)') model_dir
  read(args(3), '(a)') model_tags
  read(args(4), *) lat0
  read(args(5), *) lat1
  read(args(6), *) nlat
  read(args(7), *) lon0
  read(args(8), *) lon1
  read(args(9), *) nlon
  read(args(10), *) radius 
  read(args(11), '(a)') out_file

  !-- check validity of parameters

  !-90 <= lat0 < lat1 <=90
  if ( .not.(lat0 >= -90.d0 .and. lat0 <= 90.d0 .and. &
             lat1 >= -90.d0 .and. lat1 <= 90.d0 .and. &
             lat0 < lat1) ) then
    stop "ERROR: wrong inputs of lat0,lat1"
  endif

  !-180 <= lon0,lon1 <=180
  if ( .not.(lon0 >= -180.d0 .and. lat0 <= 180.d0 .and. &
             lon1 >= -180.d0 .and. lat1 <= 180.d0) ) then
    stop "ERROR: wrong range of lon0,lon1"
  endif
  ! in the case of great circle on longitude
  if (lon0 > lon1) then
    print *, "# lon0 > lon1, so 360.0 is added to lon1 "
    lon1 = lon1 + 360.d0
  endif

  !===== generate mesh grid of the cross-section

  ! convert units, non-dimensionalize as in SPECFEM3D_GLOBE 
  lat0 = lat0 * DEGREES_TO_RADIANS
  lat1 = lat1 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  lon1 = lon1 * DEGREES_TO_RADIANS
  radius = radius / R_EARTH_KM

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
      call geographic_geodetic2ecef(lat(ilat), lon(ilon), 0.0d0, &
                                    vr(1), vr(2), vr(3)) 
      vr = vr / sqrt(sum(vr**2))
      xyz(:,idx) = real(vr * radius, kind=CUSTOM_REAL)
    enddo
  enddo

  !===== read model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)

  !===== locate xyz in the mesh
  call sem_constants_set(iregion)

  ! initialize arrays for each mesh chunk
  allocate(model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel), &
           model_interp(nmodel,npoint), &
           uvw(3,npoint), &
           hlagrange(NGLLX,NGLLY,NGLLZ,npoint), &
           misloc(npoint), &
           misloc_final(npoint), &
           elem_ind(npoint), &
           elem_ind_final(npoint), &
           statloc(npoint), &
           statloc_final(npoint))

  ! initialize some variables
  statloc_final = -1
  misloc_final = HUGE(1.0_CUSTOM_REAL)
  elem_ind_final = -1
  model_interp = FILLVALUE

  ! typical element size at surface
  typical_size = real( max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                           ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                       * DEGREES_TO_RADIANS * R_UNIT_SPHERE, kind=CUSTOM_REAL)

  ! loop each mesh chunk
  do iproc = 0, NPROCTOT_VAL-1

    print *, "# iproc=", iproc

    call sem_mesh_read(mesh_data, mesh_dir, iproc, iregion)

    dims = [NGLLX, NGLLY, NGLLZ, mesh_data%nspec]

    call sem_io_read_gll_modeln(model_dir, iproc, iregion, &
      nmodel, model_names, dims, model_gll)

    call sem_mesh_locate_xyz(mesh_data, npoint, xyz, idoubling, &
      uvw, hlagrange, misloc, elem_ind, statloc)

    print *, "# number of located points: ", count(statloc>=0)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! safety check
      if (statloc_final(ipoint) == 1 .and. statloc(ipoint) == 1 ) then
        print *, "WARN: ipoint=", ipoint
        print *, "====: this point is located inside more than one element!"
        print *, "====: some problem may occur."
        print *, "====: only use the first located element."
        cycle
      endif

      ! for point located inside one element but not happened before 
      ! or closer to one element than located before
      if ( (statloc(ipoint) == 1 .and. statloc_final(ipoint) /= 1 ) .or. &
           (statloc(ipoint) == 0 .and. misloc(ipoint) < misloc_final(ipoint)) ) then

        ! interpolate model
        do imodel = 1, nmodel
          model_interp(imodel,ipoint) = sum( &
            hlagrange(:,:,:,ipoint) * model_gll(:,:,:,elem_ind(ipoint),imodel) )
        enddo

        statloc_final(ipoint) = statloc(ipoint)
        misloc_final(ipoint) = misloc(ipoint)
        elem_ind_final(ipoint) = elem_ind(ipoint)

      endif

    enddo ! ipoint

  enddo ! iproc

  ! convert misloc to relative to typical element size
  where (statloc_final /= -1) misloc_final = misloc_final / typical_size

  !===== output netcdf file

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

  ! define data variables: model values
  dimids = [lon_dimid, lat_dimid]
  allocate(model_varids(nmodel))
  do imodel = 1, nmodel
    call check( nf90_def_var(ncid, trim(model_names(imodel)), NF90_REAL& 
                             , dimids, model_varids(imodel)) )
    call check( nf90_put_att(ncid, model_varids(imodel) &
                             , "_FillValue", FILLVALUE) )
  enddo
  ! define data: location status and mislocation
  call check( nf90_def_var(ncid, statloc_name, NF90_INT, dimids, statloc_varid) )
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

  ! write data: interpolated model values
  allocate(model_var(nlon, nlat))
  do imodel = 1, nmodel

    ! reshape model variable into a NDIMS matrix
    do ilat = 1, nlat
      do ilon = 1, nlon
        idx = ilon + nlon*(ilat - 1)
        model_var(ilon, ilat) = model_interp(imodel, idx)
      enddo
    enddo

    call check( nf90_put_var(ncid, model_varids(imodel), model_var) )
  enddo

  ! write data: status of location and mislocation
  ! reshape variable into an NDIMS matrix
  allocate(statloc_var(nlon, nlat), misloc_var(nlon, nlat))
  do ilat = 1, nlat
    do ilon = 1, nlon
      idx = ilon + nlon*(ilat - 1)
      statloc_var(ilon, ilat) = statloc_final(idx)
      misloc_var(ilon, ilat) = misloc_final(idx)
    enddo
  enddo

  call check( nf90_put_var(ncid, statloc_varid, statloc_var) )
  call check( nf90_put_var(ncid, misloc_varid, misloc_var) )

  ! close file 
  call check( nf90_close(ncid) )


!********************
contains

subroutine check(status)
  integer, intent (in) :: status
  
  if(status /= nf90_noerr) then 
    print "('ERROR: ',a)", trim(nf90_strerror(status))
    stop "Stopped"
  endif
end subroutine check 


end program
