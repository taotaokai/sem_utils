subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_smooth_spherical_cap "
  print '(a)', "    - smooth SEM GLL model with a spherical cap kernel"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_smooth \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_tags> \"
  print '(a)', "    <sigma_h> <sigma_r> <output_dir> <out_suffix>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  The smoothing kernel is "
  print '(a)', "    ker(x,y) ~ exp(- (dtheta^2/(2*sigma_theta^2) - dr^2/(2*sigma_r^2)))"
  print '(a)', "    , where dtheta = angle between x and y, dr = radius difference between x and y "
  print '(a)', "    and sigma_theta = sigma_h/r(y)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc: number of mesh slices"
  print '(a)', "  (string) mesh_dir: mesh topology file proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir: model file proc*_reg1_<model_tags>.bin"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) sigma_h: horizontal smoothing length in km (1 std dev.)"
  print '(a)', "  (float) sigma_r: radial smoothing length in km (1 std dev.)"
  print '(a)', "  (string) output_dir: output directory of smoothed gll files"
  print '(a)', "  (string) out_suffix: proc*_reg1_<model_tags><out_suffix>.bin"
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  1) This program must run in parallel, e.g. mpirun -n <nproc> ..."

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_smooth
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  !-- command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)

  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir
  character(len=MAX_STRING_LEN) :: model_tags

  real(dp) :: sigma_h, sigma_r

  character(len=MAX_STRING_LEN) :: output_dir
  character(len=MAX_STRING_LEN) :: out_suffix

  !-- region id
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle

  !-- local variables
  integer :: i, ier, iglob, ispec
  integer :: igllx, iglly, igllz

  !-- mpi
  integer :: myrank, nrank

  !-- model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)
  character(len=MAX_STRING_LEN), allocatable :: output_model_names(:)

  !-- target mesh slice 
  ! of which the smoothed model values are calculated by averaging over the 
  ! neighbouring gll points (either in the same slice or in surrounding slices)
  type(sem_mesh_data) :: mesh_target
  integer :: iproc_target, ispec_target, nspec_target
  integer :: igllx_target, iglly_target, igllz_target
  real(dp), allocatable :: xyz_gll_target(:,:,:,:,:) ! target gll points
  real(dp), allocatable :: xyz_elem_target(:,:) ! target element centers 
  real(dp), allocatable :: model_gll_target(:,:,:,:,:) ! model gll points
  real(dp), allocatable :: weight_gll_target(:,:,:,:) ! weights on gll points

  !-- contrbuting mesh slice
  ! whose model values contribute to the target mesh model points weighted
  ! by the distance
  type(sem_mesh_data) :: mesh_contrib
  integer :: iproc_contrib, ispec_contrib, nspec_contrib
  integer :: igllx_contrib, iglly_contrib, igllz_contrib
  real(dp), allocatable :: xyz_gll_contrib(:,:,:,:,:)
  real(dp), allocatable :: xyz_elem_contrib(:,:)
  real(dp), allocatable :: volume_gll_contrib(:,:,:,:)
  real(dp), allocatable :: model_gll_contrib(:,:,:,:,:)

  !-- search parameters
  real(dp) :: sigma2_h, sigma2_r, sigma2_theta
  !real(dp) :: elem_size, max_search_dist2, dist2
  real(dp) :: max_search_dist2, dist2
  real(dp) :: dist2_radial, dist2_theta
  real(dp) :: r_target, r_contrib
  real(dp) :: xyz_target(3), xyz_contrib(3)
  ! gll points in one element
  real(dp) :: weight

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_interp_mesh: check your input arguments."
      call abort_mpi()
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
  read(args(5), *) sigma_h
  read(args(6), *) sigma_r
  read(args(7), '(a)') output_dir
  read(args(8), '(a)') out_suffix

  !===== parse model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)
  endif

  ! generate output model names
  allocate(output_model_names(nmodel))
  do imodel = 1, nmodel
    output_model_names(imodel) = trim(model_names(imodel))//trim(out_suffix)
  enddo

  !====== smoothing parameters 
  ! non-dimensionalize smoothing scales 
  ! get sigma squared
  sigma2_h = (sigma_h / EARTH_R_KM)**2
  sigma2_r = (sigma_r / EARTH_R_KM)**2

  ! typical element size at surface
  ! read from values_from_mesher.h
  !elem_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
  !                   ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
  !               * DEGREES_TO_RADIANS * R_UNIT_SPHERE

  ! maximum distance square between the contributing element center and target element center
  !max_search_dist2 = (7.0*sqrt(max(sigma2_h, sigma2_v)) + 3.0*elem_size)**2
  max_search_dist2 = 10.0*max(sigma2_h, sigma2_r)
  ! d**2 = 14*sigma**2 corresponds to exp(-14/2) = 0.0009

  !====== loop each target mesh slice
  do iproc_target = myrank, (nproc - 1), nrank

    print *, "# iproc_target=", iproc_target

    !------ read in target mesh slice
    call sem_mesh_read(mesh_dir, iproc_target, iregion, mesh_target)

    nspec_target = mesh_target%nspec

    if (allocated(xyz_elem_target)) then
      deallocate(xyz_elem_target, xyz_gll_target)
      deallocate(model_gll_target, weight_gll_target)
    endif
    allocate(xyz_elem_target(3,nspec_target))
    allocate(xyz_gll_target(3,NGLLX,NGLLY,NGLLZ,nspec_target))
    allocate(model_gll_target(nmodel,NGLLX,NGLLY,NGLLZ,nspec_target))
    allocate(weight_gll_target(NGLLX,NGLLY,NGLLZ,nspec_target))

    !------ get target xyz_gll and xyz_elem arrays
    do ispec = 1, nspec_target
      ! central points of all elements
      iglob = mesh_target%ibool(MIDX, MIDY, MIDZ, ispec)
      xyz_elem_target(:,ispec) = mesh_target%xyz_glob(:, iglob)
      ! gll points of all elements
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            iglob = mesh_target%ibool(igllx, iglly, igllz, ispec)
            xyz_gll_target(:,igllx,iglly,igllz,ispec) = mesh_target%xyz_glob(:, iglob)
          enddo
        enddo
      enddo
    enddo

    !------ collect neighbouring model values from each contributing slices
    model_gll_target = 0.0_dp
    weight_gll_target = 0.0_dp

    do iproc_contrib = 0, (nproc - 1)

      print *, "# iproc_contrib=", iproc_contrib

      !-- read contributing mesh slice
      if (iproc_contrib == iproc_target) then
        mesh_contrib = mesh_target
      else
        call sem_mesh_read(mesh_dir, iproc_contrib, iregion, mesh_contrib)
      endif

      !-- get xyz_gll and xyz_elem arrays of the contributing slice
      nspec_contrib = mesh_contrib%nspec

      if (allocated(xyz_elem_contrib)) then
        deallocate(xyz_elem_contrib, xyz_gll_contrib)
      endif
      allocate(xyz_elem_contrib(3,nspec_contrib))
      allocate(xyz_gll_contrib(3,NGLLX,NGLLY,NGLLZ,nspec_contrib))

      if (iproc_contrib == iproc_target) then
        xyz_elem_contrib = xyz_elem_target
        xyz_gll_contrib = xyz_gll_target
      else
        do ispec = 1, nspec_contrib
          ! central points of all elements 
          iglob = mesh_contrib%ibool(MIDX, MIDY, MIDZ, ispec)
          xyz_elem_contrib(:,ispec) = mesh_contrib%xyz_glob(:, iglob)
          ! gll points of all elements
          do igllz = 1, NGLLZ
            do iglly = 1, NGLLY
              do igllx = 1, NGLLX
                iglob = mesh_contrib%ibool(igllx, iglly, igllz, ispec)
                xyz_gll_contrib(:,igllx,iglly,igllz,ispec) = mesh_contrib%xyz_glob(:, iglob)
              enddo
            enddo
          enddo
        enddo
      endif

      !-- test if the distance between the contributing and target slices are more than max_search_dist2
      if (iproc_contrib /= iproc_target) then
        dist2 = huge(0.0_dp)
        do ispec_target = 1, nspec_target
          dist2 = min(dist2, minval( &
               (xyz_elem_contrib(1,:)-xyz_elem_target(1,ispec_target))**2 &
             + (xyz_elem_contrib(2,:)-xyz_elem_target(2,ispec_target))**2 &
             + (xyz_elem_contrib(3,:)-xyz_elem_target(3,ispec_target))**2))
        enddo
        if (dist2 > max_search_dist2) then
          !write(*,"(A, E12.4, A, E12.4)") "dist2= ", dist2, " is larger than max_search_dist2= ", max_search_dist2
          write(*, "(A, I4, I4, A)") "iproc_target/contrib ", iproc_target, iproc_contrib, " are far away, SKIP."
          cycle
        endif
      endif

      !-- get model_gll of the contributing slice
      if (allocated(model_gll_contrib)) then
        deallocate(model_gll_contrib)
      endif
      allocate(model_gll_contrib(nmodel, NGLLX, NGLLY, NGLLZ, nspec_contrib))
      call sem_io_read_gll_file_n(model_dir, iproc_contrib, iregion, &
        model_names, nmodel, model_gll_contrib)

      !-- calculate volume weights of each gll points in the contributing slice 
      if (allocated(volume_gll_contrib)) then
        deallocate(volume_gll_contrib)
      endif
      allocate(volume_gll_contrib(NGLLX, NGLLY, NGLLZ, nspec_contrib))
      call sem_mesh_gll_volume(mesh_contrib, volume_gll_contrib)

      !-- calculate volume weighted gaussian window average for each target element
      ! over all contributing elements
      do ispec_target = 1, nspec_target
        !write(*,"(A, I4)") "ispec_target=", ispec_target
        do ispec_contrib = 1, nspec_contrib
          !write(*,"(A, I4)") "ispec_contrib=", ispec_target

          ! calculate squared distances between the target and contrbuting elements
          dist2 = sum((xyz_elem_contrib(:,ispec_contrib)-xyz_elem_target(:,ispec_target))**2)

          ! skip if the contributing element are distant from the target element
          ! by more than max_search_dist2
          if (dist2 > max_search_dist2) then
            !write(*,"(A, E12.4, A, E12.4)") "dist2= ", dist2, " is larger than max_search_dist2= ", max_search_dist2
            cycle
          endif

          ! collect model values from the contributing element on each gll point in target element
          do igllz_target = 1, NGLLZ
            do iglly_target = 1, NGLLY
              do igllx_target = 1, NGLLX

                ! target gll point vector and radius  
                xyz_target = xyz_gll_target(:,igllx_target,iglly_target,igllz_target,ispec_target)
                r_target = sqrt(sum(xyz_target**2))

                sigma2_theta = sigma2_h/r_target**2

                do igllz_contrib = 1, NGLLZ
                  do iglly_contrib = 1, NGLLY
                    do igllx_contrib = 1, NGLLX

                      ! contributing gll point vector and radius  
                      xyz_contrib = xyz_gll_contrib(:,igllx_contrib,iglly_contrib,igllz_contrib,ispec_contrib)
                      r_contrib = sqrt(sum(xyz_contrib**2))

                      ! squared radial distance
                      dist2_radial = (r_contrib - r_target)**2
                      ! squared angular distance
                      dist2_theta = acos(sum(xyz_target*xyz_contrib)/r_target/r_contrib)

                      ! calcuate the smoothing weight on each contributing gll points
                      ! here, it is the product of smoothing weight and the gll volume (ker(x,y) * dV(y))
                      weight = exp(-0.5*dist2_radial/sigma2_r)        &
                             * exp(-0.5*dist2_theta/sigma2_theta)     &
                             * volume_gll_contrib(igllx_contrib,iglly_contrib,igllz_contrib,ispec_contrib)

                      ! add weighted model values of contributing gll point into target model gll point
                      model_gll_target(:,igllx_target,iglly_target,igllz_target,ispec_target)       &
                        = model_gll_target(:,igllx_target,iglly_target,igllz_target,ispec_target)   &
                        + weight*model_gll_contrib(:,igllx_contrib,iglly_contrib,igllz_contrib,ispec_contrib)

                      ! cumulative weight
                      weight_gll_target(igllx_target,iglly_target,igllz_target,ispec_target)  &
                        = weight_gll_target(igllx_target,iglly_target,igllz_target,ispec_target) + weight

                    enddo ! igllx_contrib = 1, NGLLX
                  enddo ! iglly_contrib = 1, NGLLY
                enddo ! igllz_contrib = 1, NGLLZ

              enddo ! igllx_target = 1, NGLLX
            enddo ! iglly_target = 1, NGLLY
          enddo ! igllz_target = 1, NGLLZ

        enddo ! ispec_contrib = 1, nspec_contrib
      enddo ! ispec_target = 1, nspec_target

    enddo ! iproc_contrib

    !-- calculate the weighted average to get the smoothed model for the target slice
    do imodel = 1, nmodel
      model_gll_target(imodel,:,:,:,:) =  model_gll_target(imodel,:,:,:,:) / weight_gll_target
    enddo

    !-- write out smoothed target model gll
    call sem_io_write_gll_file_n(output_dir, iproc_target, iregion, &
      output_model_names, nmodel, model_gll_target)

  enddo ! iproc_target

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
