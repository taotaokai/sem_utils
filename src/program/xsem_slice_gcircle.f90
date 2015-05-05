subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle - create slice of SEM model along great circle"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle <mesh_dir> <model_dir> <model_tags> "
  print '(a)', "                     <lat0> <lon0> <lat1> <lon1> <ntheta>"
  print '(a)', "                     <radius0> <radius1> <nradius> <out_file>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  create a verical slice of SEM model along specified great circle" 
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir: directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir: directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0,lon0: starting point on the great circle "
  print '(a)', "  (float) lat1,lon1: end point on the great circle "
  print '(a)', "  (integer) ntheta: number of grids along the great cicle "
  print '(a)', "  (float) radius0,radius1: begin/end radius "
  print '(a)', "  (integer) nradius: number of grids along the radius "
  print '(a)', "  (string) out_file: output file name (netcdf format) "
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_vertical_slice
  
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
  integer, parameter :: nargs = 12
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_file
  real(dp) :: lat0, lon0, lat1, lon1, radius0, radius1
  integer :: ntheta, nradius

  !-- local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier

  !-- model names
  integer :: imodel, nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- interpolation points 
  integer :: ipoint, npoint
  real(dp) :: v0(3), v1(3), v_axis(3), vr(3), rotmat(3,3)
  real(dp) :: dtheta, dradius
  real(dp), allocatable :: theta(:), radius(:)
  integer :: itheta, iradius, idx
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
  integer :: theta_dimid, radius_dimid, dimids(NDIMS)
  ! define coordinate variables/units
  integer :: theta_varid, radius_varid
  character(len=*), parameter :: UNITS = "units"
  character(len=*), parameter :: theta_name = "theta"
  character(len=*), parameter :: theta_units = "degree"
  character(len=*), parameter :: radius_name = "radius"
  character(len=*), parameter :: radius_units = "km"
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
  read(args(5), *) lon0
  read(args(6), *) lat1
  read(args(7), *) lon1
  read(args(8), *) ntheta
  read(args(9), *) radius0
  read(args(10), *) radius1
  read(args(11), *) nradius
  read(args(12), '(a)') out_file

  !===== generate mesh grid of the cross-section

  ! convert units, non-dimensionalize as in SPECFEM3D_GLOBE 
  lat0 = lat0 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  lat1 = lat1 * DEGREES_TO_RADIANS
  lon1 = lon1 * DEGREES_TO_RADIANS

  radius0 = radius0 / R_EARTH_KM
  radius1 = radius1 / R_EARTH_KM

  ! unit direction vectors at point 0,1
  call geographic_geodetic2ecef(lat0, lon0, 0.d0, v0(1), v0(2), v0(3))
  call geographic_geodetic2ecef(lat1, lon1, 0.d0, v1(1), v1(2), v1(3))

  v0 = v0 / sqrt(sum(v0**2))
  v1 = v1 / sqrt(sum(v1**2))

  ! rotation axis = v0 cross-product v1
  v_axis(1) = v0(2)*v1(3) - v0(3)*v1(2)
  v_axis(2) = v0(3)*v1(1) - v0(1)*v1(3)
  v_axis(3) = v0(1)*v1(2) - v0(2)*v1(1)
  v_axis = v_axis / sqrt(sum(v_axis**2))

  ! get grid intervals
  dtheta = acos(sum(v0*v1)) / (ntheta - 1)
  dradius = (radius1 - radius0) / (nradius - 1)

  ! mesh grid xyz(1:3, radius*theta)
  allocate(theta(ntheta), radius(nradius))
  theta = [ (i * dtheta, i=0,ntheta-1) ]
  radius = [ (radius0 + i*dradius, i=0,nradius-1 ) ]

  npoint = nradius * ntheta
  allocate(xyz(3,npoint))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    vr = matmul(rotmat, v0)
    idx = nradius*(itheta - 1)
    do iradius = 1, nradius
      xyz(1:3,idx+iradius) = REAL(radius(iradius) * vr, kind=CUSTOM_REAL)
    enddo
  enddo

  !do i = 1, npoint
  !  print *, xyz(:,i)
  !enddo

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
  model_interp = -12345.0_CUSTOM_REAL

  ! typical element size at surface
  typical_size = real( max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                           ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                       * DEGREES_TO_RADIANS * R_UNIT_SPHERE, kind=CUSTOM_REAL)

  ! loop each mesh chunk
  do iproc = 0, NPROCTOT_VAL-1

    call sem_mesh_init(mesh_data)
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_io_read_gll_modeln(model_dir, iproc, iregion, nmodel, model_names &
                                , model_gll)
    call sem_mesh_locate_xyz(mesh_data, npoint, xyz, uvw, hlagrange, misloc &
                             , elem_ind, statloc)

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

  !!debug
  !do itheta = 1, ntheta
  !  do iradius = 1, nradius
  !    idx = iradius + nradius*(itheta - 1)
  !    print *, radius(iradius)*R_EARTH_KM, theta(itheta)*RADIANS_TO_DEGREES, &
  !      model_interp(1,idx), statloc_final(idx), misloc_final(idx), &
  !      elem_ind_final(idx)
  !  enddo
  !enddo

  !===== output netcdf file
  ! create the file
  call check( nf90_create(out_file, NF90_CLOBBER, ncid))

  ! define the dimensions.
  call check( nf90_def_dim(ncid, theta_name, ntheta, theta_dimid) )
  call check( nf90_def_dim(ncid, radius_name, nradius, radius_dimid) )

  ! define the coordinate variables.
  call check( nf90_def_var(ncid, theta_name, NF90_DOUBLE, theta_dimid, theta_varid) )
  call check( nf90_def_var(ncid, radius_name, NF90_DOUBLE, radius_dimid, radius_varid) )
  ! coordinate units
  call check( nf90_put_att(ncid, theta_varid, UNITS, theta_units) )
  call check( nf90_put_att(ncid, radius_varid, UNITS, radius_units) )

  ! define data variables: model values
  dimids = [radius_dimid, theta_dimid]
  allocate(model_varids(nmodel))
  do imodel = 1, nmodel
    call check( nf90_def_var(ncid, trim(model_names(imodel)), NF90_REAL& 
                , dimids, model_varids(imodel)) )
  enddo
  ! define data: location status and mislocation
  call check( nf90_def_var(ncid, statloc_name, NF90_INT, dimids, statloc_varid) )
  call check( nf90_def_var(ncid, misloc_name, NF90_FLOAT, dimids, misloc_varid) )

  ! end define mode.
  call check( nf90_enddef(ncid) )

  ! write the coordinate variable data.
  call check( nf90_put_var(ncid, theta_varid, theta * RADIANS_TO_DEGREES) )
  call check( nf90_put_var(ncid, radius_varid, radius * R_EARTH_KM) )

  ! write data: interpolated model values
  allocate(model_var(nradius,ntheta))
  do imodel = 1, nmodel

    ! reshape model variable into a NDIMS matrix
    do itheta = 1, ntheta
      idx = nradius*(itheta - 1)
      do iradius = 1, nradius
        model_var(iradius, itheta) = model_interp(imodel, idx + iradius)
      enddo
    enddo

    call check( nf90_put_var(ncid, model_varids(imodel), model_var) )
  enddo

  ! write data: status of location and mislocation
  ! reshape variable into an NDIMS matrix
  allocate(statloc_var(nradius, ntheta), misloc_var(nradius, ntheta))
  do itheta = 1, ntheta
    idx = nradius*(itheta - 1)
    do iradius = 1, nradius
      statloc_var(iradius, itheta) = statloc_final(idx + iradius)
      misloc_var(iradius, itheta) = misloc_final(idx + iradius)
    enddo
  enddo
  call check( nf90_put_var(ncid, statloc_varid, statloc_var) )
  call check( nf90_put_var(ncid, misloc_varid, misloc_var) )

  ! close file 
  call check( nf90_close(ncid) )

!-------------------
contains

subroutine check(status)
  integer, intent (in) :: status
  
  if(status /= nf90_noerr) then 
    print "('ERROR: ',a)", trim(nf90_strerror(status))
    stop "Stopped"
  endif
end subroutine check 

end program
