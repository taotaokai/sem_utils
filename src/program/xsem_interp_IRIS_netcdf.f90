subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_interp_IRIS_netcdf "
  print '(a)', "    - interpolate IRIS netcdf model files into SEM mesh"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_interp_IRIS_netcdf \"
  print '(a)', "    <nc_file> <model_tags> "
  print '(a)', "    <mesh_dir> <nproc> "
  print '(a)', "    <flag_ellipticity> <gll_fillvalue> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) nc_file:  IRIS DMC consistent netcdf file "
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vs,vp "
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (int) flag_ellipticity: (0/1) 0: spherical; 1: WGS84 ellipsoid"
  print '(a)', "  (float) gll_fillvalue: values to use when GLL point is out of the model range"
  print '(a)', "  (string) out_dir:  output directory of new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_interp_IRIS_netcdf
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic
  use netcdf

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: nc_file, mesh_dir, model_tags, out_dir
  integer :: nproc, flag_ellipticity
  real(dp) :: gll_fillvalue 

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! model names
  integer :: nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- netcdf data
  integer :: ncid
  ! dimensions: depth,lat,lon
  integer, parameter :: NDIMS = 3
  integer :: depth_dimid, lat_dimid, lon_dimid
  ! coordinate variables
  integer :: depth_varid, lat_varid, lon_varid
  integer :: ndepth, nlat, nlon
  real(sp), allocatable :: depth_var(:), lat_var(:), lon_var(:)
  character(len=*), parameter :: depth_name = "depth"
  character(len=*), parameter :: lat_name = "latitude"
  character(len=*), parameter :: lon_name = "longitude"
  ! model variables (units is not determined)
  integer, allocatable :: model_varids(:)
  real(sp), allocatable :: model_vars(:,:,:,:)
  real(sp), allocatable :: model_fillvalues(:)

  !-- GLL mesh
  type(sem_mesh_data) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec, nspec
  ! GLL model
  real(dp), allocatable :: model_gll(:,:,:,:,:)
  ! GLL points
  real(dp) :: xyz(3)
  real(dp) :: lat, lon, depth
  ! interpolation
  integer :: imodel
  integer :: idepth_top, ilat_south, ilon_west
  real(dp) :: lat0, lat1, lon0, lon1, depth0, depth1
  real(dp) :: wx0, wx1, wy0, wy1, wz0, wz1

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_interp_IRIS_netcdf: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') nc_file
  read(args(2), '(a)') model_tags
  read(args(3), '(a)') mesh_dir
  read(args(4), *) nproc
  read(args(5), *) flag_ellipticity
  read(args(6), *) gll_fillvalue 
  read(args(7), '(a)') out_dir

  !-- validate input arguments
  if (flag_ellipticity /= 0 .and. flag_ellipticity /= 1) then
    print *, "[ERROR] flag_ellipticity must be 0/1"
    stop
  endif

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel) 

  !===== read IRIS DMC netcdf model

  call check( nf90_open(nc_file, NF90_NOWRITE, ncid) )

  ! get dimensions
  call check( nf90_inq_dimid(ncid, depth_name, depth_dimid) )
  call check( nf90_inq_dimid(ncid, lat_name, lat_dimid) )
  call check( nf90_inq_dimid(ncid, lon_name, lon_dimid) )

  call check( nf90_inquire_dimension(ncid, depth_dimid, len=ndepth) )
  call check( nf90_inquire_dimension(ncid, lat_dimid, len=nlat) )
  call check( nf90_inquire_dimension(ncid, lon_dimid, len=nlon) )

  ! get coordinate variables
  call check( nf90_inq_varid(ncid, depth_name, depth_varid) )
  call check( nf90_inq_varid(ncid, lat_name, lat_varid) )
  call check( nf90_inq_varid(ncid, lon_name, lon_varid) )

  allocate(depth_var(ndepth), lat_var(nlat), lon_var(nlon))
  call check( nf90_get_var(ncid, depth_varid, depth_var) )
  call check( nf90_get_var(ncid, lat_varid, lat_var) )
  call check( nf90_get_var(ncid, lon_varid, lon_var) )

  ! get data variables
  allocate(model_varids(nmodel), model_vars(nlon,nlat,ndepth,nmodel))
  allocate(model_fillvalues(nmodel))
  do imodel = 1, nmodel
    call check( nf90_inq_varid(ncid, model_names(imodel), model_varids(imodel)) )
    call check( nf90_get_var(ncid, model_varids(imodel), model_vars(:,:,:,imodel)) )
    call check( nf90_get_att(ncid, model_varids(imodel), "_FillValue", &
                model_fillvalues(imodel)) )
  enddo

  ! close nc file
  call check( nf90_close(ncid) )

  !===== loop each mesh slice
  do iproc = 0, (nproc - 1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! read model gll
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

    ! interploate netcdf model on each gll point
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob) * R_EARTH
            
            ! get lat,lon,depth
            if (flag_ellipticity == 1) then
              call geographic_ecef2lla(xyz(1), xyz(2), xyz(3), lat, lon, depth)
              depth = -1.0 * depth 
            else
              lon = atan2(xyz(2), xyz(1))
              lat = atan2(xyz(3), sqrt(xyz(1)**2 + xyz(2)**2) )
              depth = R_EARTH - sqrt(sum(xyz**2))
            endif

            ! change units to degrees, km
            lon = lon * RADIANS_TO_DEGREES
            lat = lat * RADIANS_TO_DEGREES
            depth = depth / 1000.0_dp

            ! get model value
            ilon_west  = maxloc(lon_var, dim=1, mask=lon_var<=lon)
            ilat_south = maxloc(lat_var, dim=1, mask=lat_var<=lat)
            idepth_top = maxloc(depth_var, dim=1, mask=depth_var<=depth)

            if (ilon_west>=1 .and. ilon_west<nlon .and. &
                ilat_south>=1 .and. ilat_south<nlat .and. &
                idepth_top>=1 .and. idepth_top<ndepth) then

              lon0 = lon_var(ilon_west)
              lon1 = lon_var(ilon_west+1)
              lat0 = lat_var(ilat_south)
              lat1 = lat_var(ilat_south+1)
              depth0 = depth_var(idepth_top)
              depth1 = depth_var(idepth_top+1)

              wx0 = (lon1 - lon)/(lon1 - lon0)
              wx1 = (lon - lon0)/(lon1 - lon0)
              wy0 = (lat1 - lat)/(lat1 - lat0)
              wy1 = (lat - lat0)/(lat1 - lat0)
              wz0 = (depth1 - depth)/(depth1 - depth0)
              wz1 = (depth - depth0)/(depth1 - depth0)

              do imodel = 1, nmodel
                if (any(model_vars(ilon_west:ilon_west+1, &
                        ilat_south:ilat_south+1,idepth_top:idepth_top+1,imodel) &
                        == model_fillvalues(imodel))) then
                  model_gll(imodel,igllx,iglly,igllz,ispec) = gll_fillvalue
                else
                  model_gll(imodel,igllx,iglly,igllz,ispec) = &
                    model_vars(ilon_west,  ilat_south,  idepth_top,  imodel)*wx0*wy0*wz0 + &
                    model_vars(ilon_west,  ilat_south,  idepth_top+1,imodel)*wx0*wy0*wz1 + &
                    model_vars(ilon_west,  ilat_south+1,idepth_top,  imodel)*wx0*wy1*wz0 + &
                    model_vars(ilon_west,  ilat_south+1,idepth_top+1,imodel)*wx0*wy1*wz1 + &
                    model_vars(ilon_west+1,ilat_south,  idepth_top,  imodel)*wx1*wy0*wz0 + &
                    model_vars(ilon_west+1,ilat_south,  idepth_top+1,imodel)*wx1*wy0*wz1 + &
                    model_vars(ilon_west+1,ilat_south+1,idepth_top,  imodel)*wx1*wy1*wz0 + &
                    model_vars(ilon_west+1,ilat_south+1,idepth_top+1,imodel)*wx1*wy1*wz1
                endif
              enddo

            else

              model_gll(:,igllx,iglly,igllz,ispec) = gll_fillvalue

            endif

          enddo
        enddo
      enddo
    enddo

    ! write out model
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, &
      model_names, nmodel, model_gll)

  enddo ! iproc

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
