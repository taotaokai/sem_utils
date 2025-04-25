subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_slab_mow_model "
  print '(a)', "    - add planar slab with metastable olivin wedge(MOW) into SEM model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_slab_model \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <slab_origin_lat> <slab_origin_lon> "
  print '(a)', "    <slab_dip> <slab_strike> "
  print '(a)', "    <slab_thick> <slab_bottom_depth> <slab_dlnv> "
  print '(a)', "    <mow_top_depth> <mow_bottom_depth> "
  print '(a)', "    <mow_xi0> <mow_xi1> <mow_xi2> "
  print '(a)', "    <mow_dlnv> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  Slab is a planar slice through the sphere, and " 
  print '(a)', "  MOW is an inverted triangular region inside the slab, and"
  print '(a)', "  forms a body of revolution " 
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) slab_origin_lat/lon:  slab origin point at surface"
  print '(a)', "  (float) slab_dip/strike/thick:  slab dip,strike angles and thickness(km)"
  print '(a)', "  (float) slab_bottom_depth:  slab bottom depth(km)"
  print '(a)', "  (float) slab_dlnv:  velocity perturbation in slab"
  print '(a)', "  (float) mow_top/bottom_depth:  top/bottom depth of MOW (km)"
  print '(a)', "  (float) mow_xi0:  distance of MOW bottom vertex below slab upper plane, "
  print '(a)', "            measured as the ratio of the slab thickness"
  print '(a)', "  (float) mw_xi1/2:  distances of MOW top two vertices"
  print '(a)', "            below slab upper plane (x1 is closer to slab upper plane), "
  print '(a)', "            measured as the ratio of the slab thickness"
  print '(a)', "  (float) mow_dlnv:  velocity perturbation in MOW"
  print '(a)', "  (string) out_dir:  output directory"
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

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 18
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_dir
  integer :: nproc
  real(dp) :: slab_origin_lat, slab_origin_lon, slab_dip, slab_strike, &
              slab_thick, slab_bottom_depth, slab_dlnv
  real(dp) :: mow_top_depth, mow_bottom_depth, &
              mow_xi0, mow_xi1, mow_xi2, mow_dlnv

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  integer :: nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  type(sem_mesh_data) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec, nspec
  real(dp), allocatable :: model_gll(:,:,:,:,:)

  real(dp) :: slab_bottom_radius, slab_upper_dist, slab_lower_dist, &
              mow_top_radius, mow_bottom_radius, &
              mow_upper_dist, mow_lower_dist, &
              theta_v0_slab_normal, &
              theta_v1_slab_normal, theta_v2_slab_normal

  real(dp), dimension(3) :: slab_origin_north, slab_origin_east, &
                            slab_origin, slab_normal, trench_normal, &
                            slab_parallel, trench_parallel, &
                            mow_v0, mow_v1, mow_v2, v01, v02, &
                            mow_upper_normal0, mow_lower_normal0, &
                            mow_upper_normal, mow_lower_normal

  real(dp) :: xyz(3), radius, ratio, v(3), vx, vy, theta, rotmat(3,3)
  real(dp) :: dist_xyz_slab_normal, dist_xyz_mow_upper_normal, &
              dist_xyz_mow_lower_normal

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_slab_mow_model: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), *) slab_origin_lat
  read(args(6), *) slab_origin_lon
  read(args(7), *) slab_dip
  read(args(8), *) slab_strike
  read(args(9), *) slab_thick
  read(args(10), *) slab_bottom_depth
  read(args(11), *) slab_dlnv
  read(args(12), *) mow_top_depth
  read(args(13), *) mow_bottom_depth
  read(args(14), *) mow_xi0
  read(args(15), *) mow_xi1
  read(args(16), *) mow_xi2
  read(args(17), *) mow_dlnv
  read(args(18), '(a)') out_dir

  !===== geometric parameters of slab and MOW
  ! use ECEF coordinate frame

  ! change unit: degree -> radians
  slab_origin_lat = slab_origin_lat * DEGREES_TO_RADIANS
  slab_origin_lon = slab_origin_lon * DEGREES_TO_RADIANS

  slab_dip = slab_dip * DEGREES_TO_RADIANS
  slab_strike = slab_strike * DEGREES_TO_RADIANS

  ! normalize length by EARTH_R_KM
  slab_thick = slab_thick / EARTH_R_KM

  ! radius of slab bottom depth
  slab_bottom_radius = 1.d0 - slab_bottom_depth / EARTH_R_KM

  print *, "# slab_bottom_radius=", slab_bottom_radius

  ! direction vector of slab origin
  call geographic_lla2ecef(slab_origin_lat, slab_origin_lon, 0.0_dp, &
                                slab_origin(1), slab_origin(2), slab_origin(3))

  slab_origin = slab_origin / sqrt(sum(slab_origin**2))

  print *, "# slab_origin=", slab_origin

  ! north and east directions at slab origin
  slab_origin_north = [ - sin(slab_origin_lat) * cos(slab_origin_lon), &
                        - sin(slab_origin_lat) * sin(slab_origin_lon), &
                          cos(slab_origin_lat) ]

  slab_origin_east =  [ - sin(slab_origin_lon), cos(slab_origin_lon), 0.0_dp ]

  ! trench normal direction at slab origin
  trench_normal = - sin(slab_strike) * slab_origin_north &
                  + cos(slab_strike) * slab_origin_east 

  ! normal vector of the slab plane
  slab_normal =   cos(slab_dip) * slab_origin &
                + sin(slab_dip) * trench_normal

  ! slab parallel direction vector
  slab_parallel =   sin(slab_dip) * slab_origin &
                  - cos(slab_dip) * trench_normal

  ! trench parallel direction vector: slab_normal x slab_parallel
  trench_parallel = [ &
    slab_normal(2)*slab_parallel(3) - slab_normal(3)*slab_parallel(2), &
    slab_normal(3)*slab_parallel(1) - slab_normal(1)*slab_parallel(3), &
    slab_normal(1)*slab_parallel(2) - slab_normal(2)*slab_parallel(1) ]

  ! distance from earth center to slab upper/lower plane
  slab_upper_dist = cos(slab_dip)
  slab_lower_dist = slab_upper_dist - slab_thick

  ! MOW top/bottom radius
  mow_top_radius = 1.0 - mow_top_depth / EARTH_R_KM
  mow_bottom_radius = 1.0 - mow_bottom_depth / EARTH_R_KM

  ! angles between MOW vertices and slab normal vector
  theta_v0_slab_normal = acos((slab_upper_dist - slab_thick*mow_xi0) &
                              / mow_bottom_radius)
  theta_v1_slab_normal = acos((slab_upper_dist - slab_thick*mow_xi1) &
                              / mow_top_radius)
  theta_v2_slab_normal = acos((slab_upper_dist - slab_thick*mow_xi2) &
                              / mow_top_radius)

  ! MOW vertices vectors 
  mow_v0 = mow_bottom_radius * &
    (  cos(slab_dip - theta_v0_slab_normal) * slab_origin &
     + sin(slab_dip - theta_v0_slab_normal) * trench_normal )
  mow_v1 = mow_top_radius * &   
    (  cos(slab_dip - theta_v1_slab_normal) * slab_origin &
     + sin(slab_dip - theta_v1_slab_normal) * trench_normal )
  mow_v2 = mow_top_radius * &
    (  cos(slab_dip - theta_v2_slab_normal) * slab_origin &
     + sin(slab_dip - theta_v2_slab_normal) * trench_normal )

  ! normal vector of MOW upper plane v0-v1 
  ! (in the plane of slab_normal and slas_origin )
  v01 = mow_v1 - mow_v0
  v01 = v01 / sqrt(sum(v01**2))

  mow_upper_normal0 = mow_v1 - sum(mow_v1*v01)*v01
  mow_upper_dist = sqrt(sum(mow_upper_normal0**2))
  mow_upper_normal0 = mow_upper_normal0 / mow_upper_dist

  ! normal vector of MOW lower plane x0-x2
  v02 = mow_v2 - mow_v0
  v02 = v02 / sqrt(sum(v02**2))

  mow_lower_normal0 = mow_v2 - sum(mow_v2*v02)*v02
  mow_lower_dist = sqrt(sum(mow_lower_normal0**2))
  mow_lower_normal0 = mow_lower_normal0 / mow_lower_dist

  !+==== parse model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel) 

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
    allocate(model_gll(nmodel, NGLLX,NGLLY,NGLLZ,nspec))

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
                                model_names, nmodel, model_gll)

    ! add slab model on each gll point 
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            ratio = 0.0_dp

            radius = sqrt(sum(xyz**2))

            if (radius >= slab_bottom_radius ) then

              ! distance projected on the normal of slab plane
              dist_xyz_slab_normal = sum(xyz * slab_normal)

              ! inside slab region
              if (dist_xyz_slab_normal <= slab_upper_dist .and. &
                  dist_xyz_slab_normal >= slab_lower_dist) & 
              then
                ratio = slab_dlnv
              endif

              ! inside MOW region
              if (radius >= mow_bottom_radius .and. radius <= mow_top_radius) &
              then

                ! compute the mow_upper/lower_normal for the current point
                v = xyz - sum(xyz*slab_normal)*slab_normal
                vx = sum(v*slab_parallel)
                vy = sum(v*trench_parallel)
                theta = atan2(vy,vx)

                call rotation_matrix(slab_normal, theta, rotmat)  

                mow_upper_normal = matmul(rotmat, mow_upper_normal0)
                mow_lower_normal = matmul(rotmat, mow_lower_normal0)

                ! distances projected on the normals of MOW upper/lower plane
                dist_xyz_mow_upper_normal = sum(xyz * mow_upper_normal)
                dist_xyz_mow_lower_normal = sum(xyz * mow_lower_normal)

                if (dist_xyz_mow_upper_normal <= mow_upper_dist .and. &
                    dist_xyz_mow_lower_normal >= mow_lower_dist) &
                then
                  ratio = mow_dlnv
                endif

              endif ! inside MOW

            endif ! above slab bottom depth

            ! get the new model
            model_gll(:, igllx, iglly, igllz, ispec) = &
              (1.0_dp + ratio) * model_gll(:, igllx, iglly, igllz, ispec)

          enddo
        enddo
      enddo
    enddo

    ! write out model
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, &
      model_names, nmodel, model_gll)

  enddo ! iproc

end program
