subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_smooth "
  print '(a)', "    - smooth SEM GLL model with a Gaussian kernel"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_smooth \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_tags> \"
  print '(a)', "    <sigma_h> <sigma_v> <output_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  Smoothing SEM model convolving with a Gaussian kernel: "
  print '(a)', "    ker(x) ~ exp(- (|x_h|^2/(2*sigma_h^2) - |x_v|^2/(2*sigma_v^2)))"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc: number of mesh slices"
  print '(a)', "  (string) mesh_dir: mesh topology file proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir: model file proc*_reg1_<model_tags>.bin"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) sigma_h: horizontal smoothing length in km (1 std dev.)"
  print '(a)', "  (float) sigma_v: vertical smoothing length in km (1 std dev.)"
  print '(a)', "  (string) output_dir: proc*_reg1_<model_tags>_smooth.bin"
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
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)

  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir
  character(len=MAX_STRING_LEN) :: model_tags

  real(dp) :: sigma_h, sigma_v 

  character(len=MAX_STRING_LEN) :: output_dir

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

  !-- target mesh slice 
  ! of which the smoothed model values are calculated by averaging over the 
  ! neighbouring gll points (either in the same slice or in surrounding slices)
  type(sem_mesh_data) :: mesh_target
  integer :: iproc_target, ispec_target, nspec_target
  real(dp), allocatable :: xyz_gll_target(:,:,:,:,:) ! target gll points
  real(dp), allocatable :: xyz_elem_target(:,:) ! target element centers 
  real(dp), allocatable :: model_gll_target(:,:,:,:,:) ! model gll points
  real(dp), allocatable :: weight_gll_target(:,:,:,:) ! weight gll points

  !-- contrbuting mesh slice
  ! whose model values contribute to the target mesh model points weighted
  ! by the distance
  type(sem_mesh_data) :: mesh_contrib
  integer :: iproc_contrib, ispec_contrib, nspec_contrib
  real(dp), allocatable :: xyz_contrib(:,:,:,:,:) ! contributing points
  real(dp), allocatable :: xyz_gll_contrib(:,:,:,:,:)
  real(dp), allocatable :: xyz_elem_contrib(:,:)
  real(dp), allocatable :: volume_gll_contrib(:,:,:,:)
  real(dp), allocatable :: model_gll_contrib(:,:,:,:,:)
  real(dp), allocatable :: weight_gll_contrib(:,:,:,:)

  !-- search parameters
  real(dp) :: sigma2_h, sigma2_v
  !real(dp) :: elem_size, max_search_dist2, dist2
  real(dp) :: max_search_dist2, dist2
  real(dp) :: r_unit_target(3)
  real(dp), dimension(3,NGLLX,NGLLY,NGLLZ) :: xyz_gll
  real(dp), dimension(NGLLX,NGLLY,NGLLZ) :: dist2_v_gll, dist2_h_gll, gauss_weight_gll

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
  read(args(5), '(a)') sigma_h
  read(args(6), '(a)') sigma_v
  read(args(7), '(a)') output_dir

  !===== parse model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)
  endif

  !====== smoothing parameters 
  ! non-dimensionalize smoothing scales 
  ! get sigma squared
  sigma2_h = (sigma_h / R_EARTH_KM)**2
  sigma2_v = (sigma_v / R_EARTH_KM)**2

  ! typical element size at surface
  ! read from values_from_mesher.h
  !elem_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
  !                   ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
  !               * DEGREES_TO_RADIANS * R_UNIT_SPHERE

  ! maximum distance square between the contributing element center and target element center
  !max_search_dist2 = (7.0*sqrt(max(sigma2_h, sigma2_v)) + 3.0*elem_size)**2
  max_search_dist2 = 10.0*max(sigma2_h, sigma2_v)
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
    allocate(model_gll_target(3,NGLLX,NGLLY,NGLLZ,nspec_target))
    allocate(weight_gll_target(NGLLX,NGLLY,NGLLZ,nspec_target))

    !------ get target xyz_gll and xyz_elem arrays
    do ispec = 1, nspec_target
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            ! central points of all elements 
            iglob = mesh_target%ibool(MIDX, MIDY, MIDZ, ispec)
            xyz_elem_target(:,ispec) = mesh_target%xyz_glob(:, iglob)
            ! gll points of all elements
            iglob = mesh_target%ibool(igllx, iglly, igllz, ispec)
            xyz_gll_target(:,igllx,iglly,igllz,ispec) = mesh_target%xyz_glob(:, iglob)
          enddo
        enddo
      enddo
    enddo

    !------ loop each contributing slices
    model_gll_target = 0.0_dp
    weight_gll_target = 0.0_dp

    do iproc_contrib = 0, (nproc - 1)

      print *, "# iproc_contrib=", iproc_contrib

      !-- read contributing mesh slice
      call sem_mesh_read(mesh_dir, iproc_contrib, iregion, mesh_contrib)

      !-- get xyz_gll and xyz_elem arrays of the contributing slice
      nspec_contrib = mesh_contrib%nspec

      if (allocated(xyz_elem_contrib)) then
        deallocate(xyz_elem_contrib, xyz_gll_contrib)
      endif
      allocate(xyz_elem_contrib(3,nspec_contrib))
      allocate(xyz_gll_contrib(3,NGLLX,NGLLY,NGLLZ,nspec_contrib))

      do ispec = 1, nspec_contrib
        do igllz = 1, NGLLZ
          do iglly = 1, NGLLY
            do igllx = 1, NGLLX
              ! central points of all elements 
              iglob = mesh_contrib%ibool(MIDX, MIDY, MIDZ, ispec)
              xyz_elem_contrib(:,ispec) = mesh_contrib%xyz_glob(:, iglob)
              ! gll points of all elements
              iglob = mesh_contrib%ibool(igllx, iglly, igllz, ispec)
              xyz_gll_contrib(:,igllx,iglly,igllz,ispec) = mesh_contrib%xyz_glob(:, iglob)
            enddo
          enddo
        enddo
      enddo

      !-- test if the distance between the contributing and target slices are more than 
      ! max_search_dist2
      dist2 = huge(0.0_dp)
      do ispec = 1, nspec_target
        dist2 = min(dist2, &
          minval( (xyz_elem_contrib(1,:) - xyz_elem_target(1,ispec))**2 &
                + (xyz_elem_contrib(2,:) - xyz_elem_target(2,ispec))**2 &
                + (xyz_elem_contrib(3,:) - xyz_elem_target(3,ispec))**2))
      enddo
      if (dist2 > max_search_dist2) then
        print *, "# [INFO] contributing and target slices are too far away, skip."
        cycle
      endif

      !-- get model_gll of the contributing slice
      if (allocated(model_gll_contrib)) then
        deallocate(model_gll_contrib)
      endif
      allocate(model_gll_contrib(nmodel, NGLLX, NGLLY, NGLLZ, nspec_contrib))
      call sem_io_read_gll_file_n(model_dir, iproc_contrib, iregion, &
        model_names, nmodel, model_gll_contrib)

      !-- calculate volume weights of each gll points for the contributing slice 
      if (allocated(volume_gll_contrib)) then
        deallocate(volume_gll_contrib)
      endif
      allocate(volume_gll_contrib(NGLLX, NGLLY, NGLLZ, nspec_contrib))
      call sem_mesh_gll_volume(mesh_contrib, volume_gll_contrib)

      !-- calculate volume weighted gaussian window average for each target element
      ! over all contributing elements
      do ispec_target = 1, nspec_target
        do ispec_contrib = 1, nspec_contrib

          ! calculate squared distances between the target and contrbuting elements
          dist2 = sum(xyz_elem_contrib(:,ispec_contrib) - xyz_elem_target(:,ispec_target)**2)

          ! skip if the contributing element are distant from the target element
          ! by more than max_search_dist2
          if (dist2 > max_search_dist2) then
            cycle
          endif

          ! loop each gll point in target element
          do igllz = 1, NGLLZ
            do iglly = 1, NGLLY
              do igllx = 1, NGLLX
                ! unit vertical(radial) vector through the target gll point  
                r_unit_target = xyz_gll_target(:,igllx,iglly,igllz,ispec_target)
                r_unit_target = r_unit_target/sqrt(sum(r_unit_target**2))
                ! horizontal and vertical distance squared from contributing gll points to the target gll point
                xyz_gll(1,:,:,:) = xyz_gll_contrib(1,:,:,:,ispec_contrib) - xyz_gll_target(1,igllx,iglly,igllz,ispec_target)
                xyz_gll(2,:,:,:) = xyz_gll_contrib(2,:,:,:,ispec_contrib) - xyz_gll_target(2,igllx,iglly,igllz,ispec_target)
                xyz_gll(3,:,:,:) = xyz_gll_contrib(3,:,:,:,ispec_contrib) - xyz_gll_target(3,igllx,iglly,igllz,ispec_target)
                dist2_v_gll = (  xyz_gll(1,:,:,:)*r_unit_target(1) &
                               + xyz_gll(2,:,:,:)*r_unit_target(2) &
                               + xyz_gll(3,:,:,:)*r_unit_target(3) )**2
                dist2_h_gll = sum(xyz_gll**2, dim=1) - dist2_v_gll
                ! calcuate the smoothing weight on each contributing gll points
                gauss_weight_gll = exp(-0.5*dist2_v_gll/sigma2_v - 0.5*dist2_h_gll/sigma2_h)
                ! add to target model gll point
                do imodel = 1, nmodel
                  model_gll_target(imodel,igllx,iglly,igllz,ispec_target) = &
                    model_gll_target(imodel,igllx,iglly,igllz,ispec_target) &
                    + sum(model_gll_contrib(imodel,:,:,:,ispec_contrib) &
                         * volume_gll_contrib(:,:,:,ispec_contrib) &
                         * gauss_weight_gll)
                  weight_gll_target(igllx,iglly,igllz,ispec_target) = &
                    weight_gll_target(igllx,iglly,igllz,ispec_target) &
                    + sum(volume_gll_contrib(:,:,:,ispec_contrib) * gauss_weight_gll)
                enddo
              enddo
            enddo
          enddo

        enddo !do ispec_contrib = 1, nspec_contrib
      enddo !do ispec_target = 1, nspec_target

    enddo ! iproc_contrib

    !-- get the smoothed model for the target slice
    do imodel = 1, nmodel
      model_gll_target(imodel,:,:,:,:) =  model_gll_target(imodel,:,:,:,:) / weight_gll_target
    enddo

    !-- write out smoothed target model gll
    do imodel = 1, nmodel
      model_names(imodel) = trim(model_names(imodel))//'_smooth'
    enddo
    call sem_io_write_gll_file_n(output_dir, iproc_target, iregion, &
      model_names, nmodel, model_gll_target)

  enddo ! iproc_target

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
