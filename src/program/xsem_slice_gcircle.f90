subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle - create slice of SEM model along great circle"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <lat0> <lon0> <azimuth> "
  print '(a)', "    <theta0> <theta1> <ntheta> "
  print '(a)', "    <radius0> <radius1> <nradius> "
  print '(a)', "    <out_file> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  create a verical slice of SEM model along specified great circle" 
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir: directory containing proc*_reg1_solver_data.bin"
  print '(a)', "  (int) nproc: number of mesh slices"
  print '(a)', "  (string) model_dir: directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0,lon0: origin point on the great circle (degrees)"
  print '(a)', "  (float) azimuth: shooting azimuth of the great circle (degrees)"
  print '(a)', "  (float) theta0, theta1: begin/end angles along great circle (degrees) "
  print '(a)', "  (integer) ntheta: number of grids along the great cicle "
  print '(a)', "  (float) radius0,radius1: begin/end radius (km) "
  print '(a)', "  (integer) nradius: number of grids along the radius "
  print '(a)', "  (string) out_file: output file name (GMT netcdf format) "
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

  !-- command line args
  integer, parameter :: nargs = 14
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_file
  real(dp) :: lat0, lon0, azimuth, theta0, theta1, radius0, radius1
  integer :: nproc, ntheta, nradius

  !-- local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  real(sp), parameter :: FILLVALUE_sp = huge(0._sp)
  integer :: i, iproc, ier

  !-- model names
  integer :: imodel, nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- interpolation points 
  integer :: ipoint, npoint
  real(dp) :: v0(3), v1(3), vnorth(3), veast(3), v_axis(3), vr(3), rotmat(3,3)
  real(dp) :: dtheta, dradius
  real(dp), allocatable :: theta(:), radius(:)
  integer :: itheta, iradius, idx
  real(dp), allocatable :: xyz(:,:)
  integer, allocatable :: idoubling(:)

  !-- sem location 
  type(sem_mesh_data) :: mesh_data
  type(sem_mesh_location), allocatable :: location_1slice(:)
  integer, parameter :: nnearest = 10
  real(dp) :: typical_size, max_search_dist, max_misloc
  real(dp), allocatable :: final_misloc(:)
  integer, allocatable :: final_stat(:)
  !-- model interpolation
  integer :: nspec
  real(dp), allocatable :: model_gll(:,:,:,:,:), model_interp(:,:)

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
  real(sp), allocatable :: model_var(:,:)
  ! define data: location status, mislocation
  integer :: stat_varid, misloc_varid
  integer, allocatable :: stat_var(:,:)
  real(sp), allocatable :: misloc_var(:,:)
  character(len=*), parameter :: stat_name = "location_status"
  character(len=*), parameter :: misloc_name = "misloc"

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    stop "[ERROR] xsem_slice_gcircle: check your input arguments."
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
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
  read(args(14), '(a)') out_file

  !===== generate mesh grid of the cross-section

  ! convert units, non-dimensionalize as in SPECFEM3D_GLOBE 
  lat0 = lat0 * DEGREES_TO_RADIANS
  lon0 = lon0 * DEGREES_TO_RADIANS
  azimuth = azimuth * DEGREES_TO_RADIANS

  theta0 = theta0 * DEGREES_TO_RADIANS
  theta1 = theta1 * DEGREES_TO_RADIANS

  radius0 = radius0 / R_EARTH_KM
  radius1 = radius1 / R_EARTH_KM

  ! unit direction vector v0 at origin point of the great circle
  call geographic_geodetic2ecef(lat0, lon0, 0.d0, v0(1), v0(2), v0(3))
  v0 = v0 / sqrt(sum(v0**2))

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
  allocate(theta(ntheta), radius(nradius))
  theta = [ (theta0 + i * dtheta, i=0,ntheta-1) ]
  radius = [ (radius0 + i * dradius, i=0,nradius-1 ) ]

  npoint = nradius * ntheta
  allocate(xyz(3, npoint))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    vr = matmul(rotmat, v0)
    idx = nradius * (itheta - 1)
    do iradius = 1, nradius
      xyz(:, idx+iradius) = radius(iradius) * vr
    enddo
  enddo

  ! layer ID 
  allocate(idoubling(npoint))
  idoubling = IFLAG_DUMMY

  !===== parse model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)

  !===== locate xyz in each mesh slice

  !-- initialize variables
  allocate(location_1slice(npoint))

  allocate(final_stat(npoint), final_misloc(npoint))
  final_stat = -1
  final_misloc = HUGE(1.0_dp)

  allocate(model_interp(nmodel, npoint))
  model_interp = FILLVALUE_sp

  ! typical element size at surface
  typical_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                     ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                 * DEGREES_TO_RADIANS * R_UNIT_SPHERE

  max_search_dist = 5.0 * typical_size

  max_misloc = typical_size / 4.0

  !-- loop each mesh slice
  loop_proc: do iproc = 0, (nproc - 1)

    print *, "# iproc=", iproc

    ! read mesh
    call sem_mesh_read(mesh_data, mesh_dir, iproc, iregion)

    ! read model
    nspec = mesh_data%nspec
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(nmodel, NGLLX, NGLLY, NGLLZ, nspec))

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
      nmodel, model_names, model_gll)

    ! locate points in this mesh slice
    call sem_mesh_locate_kdtree2(mesh_data, npoint, xyz, idoubling, &
      nnearest, max_search_dist, max_misloc, location_1slice)

    print *, "# number of located points: ", count(location_1slice%stat>=0)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! safety check
      if (final_stat(ipoint) == 1 .and. &
          location_1slice(ipoint)%stat == 1 ) then
        print *, "[WARN] ipoint=", ipoint
        print *, "------ this point is located inside more than one element!"
        print *, "------ some problem may occur."
        print *, "------ only use the first located element."
        cycle
      endif

      ! for point located inside one element but not happened before 
      ! or closer to one element than located before
      if (location_1slice(ipoint)%stat == 1 &
          .or. &
          (location_1slice(ipoint)%stat == 0 .and. &
           location_1slice(ipoint)%misloc < final_misloc(ipoint)) ) &
      then

        ! interpolate model
        do imodel = 1, nmodel
          model_interp(imodel, ipoint) = &
            sum( location_1slice(ipoint)%lagrange * &
                 model_gll(imodel, :, :, :, location_1slice(ipoint)%eid) )
        enddo

        final_stat(ipoint)   = location_1slice(ipoint)%stat
        final_misloc(ipoint) = location_1slice(ipoint)%misloc

      endif

    enddo ! ipoint

  enddo loop_proc 

  ! convert misloc to relative to typical element size
  where (final_stat /= -1) final_misloc = final_misloc / typical_size

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
  dimids = [theta_dimid, radius_dimid]
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
  call check( nf90_put_var(ncid, theta_varid, theta * RADIANS_TO_DEGREES) )
  call check( nf90_put_var(ncid, radius_varid, radius * R_EARTH_KM) )

  ! write data: interpolated model values
  allocate(model_var(ntheta, nradius))
  do imodel = 1, nmodel

    ! reshape model variable into a NDIMS matrix
    do itheta = 1, ntheta
      idx = nradius * (itheta - 1)
      do iradius = 1, nradius
        model_var(itheta, iradius) = &
          real(model_interp(imodel, idx + iradius), kind=sp)
      enddo
    enddo

    call check( nf90_put_var(ncid, model_varids(imodel), model_var) )
  enddo

  ! write data: status of location and mislocation
  ! reshape variable into an NDIMS matrix
  allocate(stat_var(ntheta, nradius), misloc_var(ntheta, nradius))
  do itheta = 1, ntheta
    idx = nradius*(itheta - 1)
    do iradius = 1, nradius
      stat_var(itheta, iradius) = final_stat(idx + iradius)
      misloc_var(itheta, iradius) = real(final_misloc(idx + iradius), kind=sp)
    enddo
  enddo

  call check( nf90_put_var(ncid, stat_varid, stat_var) )
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
