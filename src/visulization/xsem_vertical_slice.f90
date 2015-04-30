subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle - create slice of SEM model along great circle"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slice_gcircle <mesh_dir> <model_dir> <model_tags> "
  print '(a)', "                     <lat0> <lon0> <lat1> <lon1> <ntheta>"
  print '(a)', "                     <r0> <r1> <nr> <out_file>"
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
  print '(a)', "  (float) r0,r1: start/end radius "
  print '(a)', "  (integer) nr: number of grids along the radius "
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

  !---- declare variables
  ! command line args
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, xyz_list, out_list

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier 

  integer :: imodel, nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  integer :: ipoint, npoint
  real(CUSTOM_REAL), allocatable :: xyz(:,:)

  type (mesh) :: mesh_data
  integer :: igllx, iglly, igllz
  real(CUSTOM_REAL), allocatable :: model_gll(:,:,:,:,:), model_interp(:,:)
  real(CUSTOM_REAL), allocatable :: uvw(:,:), hlagrange(:,:,:,:), misloc(:)
  integer, allocatable :: elem_ind(:), stat_loc(:), stat_loc_final(:)

  !---- read command line arguments
  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (i <= nargs .and. len_trim(args(i)) == 0) then
      call selfdoc()
      stop
    endif
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), '(a)') model_dir
  read(args(3), '(a)') model_tags
  read(args(4), '(a)') xyz_list
  read(args(5), '(a)') out_list

  !---- generate mesh grid of the cross-section
  call geodetic2ecef(lat0,lon0,0.0,x0,y0,z0)
  call geodetic2ecef(lat1,lon1,0.0,x1,y1,z1)

  ! normalize v0,v1 to uint vector 
  norm_v0 = sqrt(sum(x0**2 + y0**2 + z0**2))
  v0(1) = x0 / norm_v0
  v0(2) = y0 / norm_v0
  v0(3) = z0 / norm_v0

  norm_v1 = sqrt(sum(x1**2 + y1**2 + z1**2))
  v1(1) = x1 / norm_v1
  v1(2) = y1 / norm_v1
  v1(3) = z1 / norm_v1

  ! rotation axis = v0 cross-product v1
  v_axis(1) = v0(2)*v1(3) - v0(3)*v1(2)
  v_axis(2) = v0(3)*v1(1) - v0(1)*v1(3)
  v_axis(3) = v0(1)*v1(2) - v0(2)*v1(1)

  ! get gird intervals
  dtheta = acos(dot_product(v0,v1)) / (ntheta - 1)
  dradius = (r0 - r1) / (nradius - 1)

  ! get mesh grid
  theta = (/(i * detha, i=0,nthta-1 )/)
  radius = (/(r0 + i*dradius, i=0,nradius-1 )/)

  npoint = nradius * ntheta
  allocate(xyz(3,npoint))

  do itheta = 1, ntheta
    call rotation_matrix(v_axis, theta(itheta), rotmat)
    do iradius = 1, nradius
      idx = iradius + nradius*(itheta - 1)
      xyz(:,idx) = radius(iradius) * matmul(rotmat, v0)
    enddo
  enddo

  !---- read model tags 
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print '(a)', '# nmodel=', nmodel
  print '(a)', '# model_names=', model_names

  !---- locate xyz in the mesh 
  call sem_constants_set(iregion)

  ! initialize arrays for each mesh chunk
  allocate(model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel), &
           model_interp(nmodel,npoint), &
           uvw(NDIM,npoint), &
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

  !---- output netcdf file

end program
