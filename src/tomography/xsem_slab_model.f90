subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slab_model - add slab/metastable olivin wedge into SEM model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slab_model <mesh_dir> <model_dir> <model_tags> "
  print '(a)', "                  <lat0> <lon0> <dip> <strike> <thickness> "
  print '(a)', "                  <dlnv_slab> <mow_w0> <mow_w1> <dlnv_mow> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  the slab model is planar, invariant along strike" 
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) lat0,lon0: trench reference point"
  print '(a)', "  (float) dip, strike, thickness: slab dip,strike and thickness"
  print '(a)', "  (float) dlnv_slab: slab velocity perturbation"
  print '(a)', "  (float) mow_w0/w1: ratios of wedge width at 410km, bottom vertex at 660km"
  print '(a)', "  (float) dlnv_mow: velocity perturbation of metastable olivine wedge"
  print '(a)', "  (string) out_dir: output directory name"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_vertical_slice
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !---- declare variables
  ! command line args
  integer, parameter :: nargs = 13
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_dir
  double precision :: lat0, lon0, height0, dip, strike, thickness, dlnv_slab 
  double precision :: mow_w0, mow_w1, dlnv_mow

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier

  integer :: nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  type (mesh) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec
  real(CUSTOM_REAL), allocatable :: model_gll(:,:,:,:,:)

  double precision :: x, y, z, north, east, down, r, d, radius, ratio
  double precision :: mow_d0, mow_d1, mow_radius0, mow_radius1, mow_dlen
  double precision :: slab_rlen, mow_r0_upper, mow_r0_lower, mow_r1, mow_r0 

  !---- read command line arguments
  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (len_trim(args(i)) == 0) then
      call selfdoc()
      stop
    endif
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), '(a)') model_dir
  read(args(3), '(a)') model_tags
  read(args(4), *) lat0
  read(args(5), *) lon0
  read(args(6), *) dip
  read(args(7), *) strike
  read(args(8), *) thickness
  read(args(9), *) dlnv_slab
  read(args(10), *) mow_w0
  read(args(11), *) mow_w1
  read(args(12), *) dlnv_mow
  read(args(13), '(a)') out_dir

  dip = dip * DEGREES_TO_RADIANS
  strike = strike * DEGREES_TO_RADIANS
  thickness = thickness * 1000.d0
  slab_rlen = thickness / sin(dip)

  mow_d0 = 410.d0 * 1000.d0 
  mow_d1 = 660.d0 * 1000.d0
  mow_radius0 = R_EARTH - mow_d1
  mow_radius1 = R_EARTH - mow_d0
  mow_dlen = mow_d1 - mow_d0
  mow_r0_upper = mow_d0 / tan(dip)
  mow_r0_lower = mow_r0_upper - slab_rlen*mow_w0
  mow_r1 = mow_d1/tan(dip) - slab_rlen*mow_w1

  !---- read model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print '(a)', '# nmodel=', nmodel
  print '(a)', '# model_names=', model_names

  !---- locate xyz in the mesh 
  call sem_constants_set(iregion)

  ! initialize arrays for each mesh chunk
  allocate(model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel))

  ! loop each mesh chunk
  do iproc = 0, NPROCTOT_VAL-1

    print '(a)', '# iproc=', iproc

    call sem_mesh_init(mesh_data)
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_io_read_gll_modeln(model_dir, iproc, iregion, nmodel, model_names &
                                , model_gll)

    ! add slab model on each gll point 
    do ispec = 1, NSPEC
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            x = mesh_data%xyz(1,iglob) * R_EARTH
            y = mesh_data%xyz(2,iglob) * R_EARTH 
            z = mesh_data%xyz(3,iglob) * R_EARTH
            radius = sqrt(x**2 + y**2 + z**2)

            ! calculate the local north,east,down coordinates relative to trench
            ! reference point
            call geographic_ecef2ned(x,y,z,lat0,lon0,height0,north,east,down)

            ! get the distance from trench/slab
            r = cos(strike)*east - sin(strike)*north
            d = (down - tan(dip)*r) / cos(dip)

            ratio = 0.0d0

            ! assign velocity perturbation of slab
            if (radius > mow_radius0 .and. d >= 0.d0 .and. d <= thickness) then
              ratio = dlnv_slab
            endif

            ! assign velocity perturbation for metastable olive wedge
            mow_r1 = mow_r0_upper + (d - mow_d0)/mow_dlen * (mow_r1 - mow_r0_upper)
            mow_r0 = mow_r0_lower + (d - mow_d0)/mow_dlen * (mow_r1 - mow_r0_lower)

            if (radius > mow_radius0 .and. radius < mow_radius1 &
                .and. r >= mow_r0 .and. r <= mow_r1) then
              ratio = dlnv_mow
            endif

            ! get the new model
            model_gll(igllx,iglly,igllz,ispec,:) = REAL((1.0d0 + ratio) * &
              model_gll(igllx,iglly,igllz,ispec,:), kind=CUSTOM_REAL)

          enddo
        enddo
      enddo
    enddo

    ! write out model
    call sem_io_write_gll_modeln(out_dir, iproc, iregion, nmodel, model_names &
                                 , model_gll)

  enddo ! iproc


end program
