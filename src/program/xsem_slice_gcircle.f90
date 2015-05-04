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
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
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
  real(dp) :: x0, y0, z0, x1, y1, z1
  real(dp) :: v0(3), v1(3), v_axis(3), rotmat(3,3)
  real(dp) :: norm_v0, norm_v1
  real(dp) :: dtheta, dradius
  real(dp), allocatable :: theta(:), radius(:)
  integer :: itheta, iradius, idx
  real(CUSTOM_REAL), allocatable :: xyz(:,:)

  !-- sem interpolation
  type (mesh) :: mesh_data
  integer :: igllx, iglly, igllz
  real(CUSTOM_REAL), allocatable :: model_gll(:,:,:,:,:), model_interp(:,:)
  real(CUSTOM_REAL), allocatable :: uvw(:,:), hlagrange(:,:,:,:), misloc(:)
  integer, allocatable :: elem_ind(:), stat_loc(:), stat_loc_final(:)

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
  ! data units is again not determined
  integer, allocatable :: model_varids(:)
  real(CUSTOM_REAL), allocatable :: model_var(:,:)

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
  radius0 = radius0 / R_EARTH_KM
  radius1 = radius1 / R_EARTH_KM

  ! unit direction vectors at point 0,1
  call geographic_geodetic2ecef(lat0, lon0, 0.d0, x0, y0, z0)
  call geographic_geodetic2ecef(lat1, lon1, 0.d0, x1, y1, z1)

  norm_v0 = sqrt(x0**2 + y0**2 + z0**2)
  v0(1) = x0 / norm_v0
  v0(2) = y0 / norm_v0
  v0(3) = z0 / norm_v0

  norm_v1 = sqrt(x1**2 + y1**2 + z1**2)
  v1(1) = x1 / norm_v1
  v1(2) = y1 / norm_v1
  v1(3) = z1 / norm_v1

  ! rotation axis = v0 cross-product v1
  v_axis(1) = v0(2)*v1(3) - v0(3)*v1(2)
  v_axis(2) = v0(3)*v1(1) - v0(1)*v1(3)
  v_axis(3) = v0(1)*v1(2) - v0(2)*v1(1)

  ! get grid intervals
  dtheta = acos(dot_product(v0,v1)) / (ntheta - 1)
  dradius = (radius0 - radius1) / (nradius - 1)

  ! mesh grid xyz(:, radius+theta)
  allocate(theta(ntheta), radius(nradius))
  theta = (/(i * dtheta, i=0,ntheta-1 )/)
  radius = (/(radius0 + i*dradius, i=0,nradius-1 )/)

  npoint = nradius * ntheta
  allocate(xyz(3,npoint))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    do iradius = 1, nradius
      idx = iradius + nradius*(itheta - 1)
      xyz(:,idx) = radius(iradius) * matmul(rotmat, v0)
    enddo
  enddo

  !===== read model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print '(a)', '# nmodel=', nmodel
  print '(a)', '# model_names=', model_names

  !===== locate xyz in the mesh
  call sem_constants_set(iregion)

  ! initialize arrays for each mesh chunk
  allocate(model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel), &
           model_interp(nmodel,npoint), &
           uvw(3,npoint), &
           hlagrange(NGLLX,NGLLY,NGLLZ,npoint), &
           misloc(npoint), &
           elem_ind(npoint), &
           stat_loc(npoint), &
           stat_loc_final(npoint))

  stat_loc_final = -1
  model_interp = 0.0_CUSTOM_REAL

  ! loop each mesh chunk
  do iproc = 0, NPROCTOT_VAL-1
  !do iproc = 0, 0

    print '(a)', '# iproc=', iproc

    call sem_mesh_init(mesh_data)
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_io_read_gll_modeln(model_dir, iproc, iregion, nmodel, model_names &
                                , model_gll)
    call sem_mesh_locate_xyz(mesh_data, npoint, xyz, uvw, hlagrange, misloc &
                             , elem_ind, stat_loc)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! for point located inside an element
      if (stat_loc(ipoint) == 1) then

        ! should not happen unless exactly on the element surface
        if (stat_loc_final(ipoint) == 1) then
          print '(a)', "WARN: more than one elements locate point ", ipoint
          print '(a)', "WARN: only use the first located element."
          cycle
        endif

        ! interpolate model
        do imodel = 1, nmodel
          model_interp(imodel,ipoint) = sum( &
            hlagrange(:,:,:,ipoint) * model_gll(:,:,:,elem_ind(ipoint),imodel) )
        enddo

        stat_loc_final(ipoint) = 1

      endif

    enddo ! ipoint

  enddo ! iproc

  !===== output netcdf file
  ! create the file
  call check( nf90_create(out_file, nf90_noclobber, ncid))

  ! define the dimensions.
  call check( nf90_def_dim(ncid, theta_name, ntheta, theta_dimid) )
  call check( nf90_def_dim(ncid, radius_name, nradius, radius_dimid) )

  ! define the coordinate variables.
  call check( nf90_def_var(ncid, theta_name, NF90_DOUBLE, theta_dimid, theta_varid) )
  call check( nf90_def_var(ncid, radius_name, NF90_DOUBLE, radius_dimid, radius_varid) )
  ! coordinate units
  call check( nf90_put_att(ncid, theta_varid, UNITS, theta_units) )
  call check( nf90_put_att(ncid, radius_varid, UNITS, radius_units) )

  ! define data variables
  dimids = [radius_dimid, theta_dimid]
  allocate(model_varids(nmodel))
  do imodel = 1, nmodel
    call check( nf90_def_var(ncid, trim(model_names(imodel)), NF90_REAL& 
                , dimids, model_varids(imodel)) )
  enddo

  ! end define mode.
  call check( nf90_enddef(ncid) )

  ! write the coordinate variable data.
  call check( nf90_put_var(ncid, theta_varid, theta) )
  call check( nf90_put_var(ncid, radius_varid, radius) )

  ! write data varaibles
  allocate(model_var(nradius,ntheta))
  do imodel = 1, nmodel

    ! reshape model variable into a NDIMS matrix
    do itheta = 1, ntheta
      do iradius = 1, nradius
        idx = iradius + nradius*(itheta - 1)
        model_var(iradius, itheta) = model_interp(imodel, idx)
      enddo
    enddo

    call check( nf90_put_var(ncid, model_varids(imodel), model_var) )

  enddo

  ! close file 
  call check( nf90_close(ncid) )

!-------------------
contains

subroutine check(status)
  integer, intent (in) :: status
  
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  endif
end subroutine check 

end program
