subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_interp_mesh "
  print '(a)', "    - interpolate GLL model from one SEM mesh onto a new mesh"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_interp_mesh \ "
  print '(a)', "    <old_mesh_dir> <nproc_old> <old_model_dir> <model_tags> "
  print '(a)', "    <new_mesh_dir> <nproc_new> <output_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  interpolate GLL model given on one SEM mesh onto  another SEM mesh"
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) old_mesh_dir: directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (int) nproc_old: number of slices of the old mesh"
  print '(a)', "  (string) old_model_dir: directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (string) new_mesh_dir: directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (int) nproc_new: number of slices of the new mesh"
  print '(a)', "  (string) output_dir: output directory for interpolated model files"
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  1) This program must run in parallel, e.g. mpirun -n <nproc> ..."
  print '(a)', "  2) use values_from_mesher.h created for the old mesh to compile"
  print '(a)', "  3) the new_model is used as the background model for the place "
  print '(a)', "     old mesh doesn't include"
! print '(a)', "  2) max_misloc ~ 0.25 * typical element size of old mesh "
! print '(a)', "  3) max_search_dist ~ 5.0 * typical element size of old mesh"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_interp_mesh
  
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
  character(len=MAX_STRING_LEN) :: old_mesh_dir, old_model_dir
  integer :: nproc_old
  character(len=MAX_STRING_LEN) :: model_tags

  character(len=MAX_STRING_LEN) :: new_mesh_dir
  integer :: nproc_new

  character(len=MAX_STRING_LEN) :: output_dir

  !-- region id
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle

  !-- local variables
  integer :: i, ier, iglob, ispec

  !-- mpi
  integer :: myrank, nrank

  !-- model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- old mesh slice
  type(sem_mesh_data) :: mesh_old
  integer :: iproc_old, nspec_old
  real(dp), allocatable :: model_gll_old(:,:,:,:,:)

  !-- new mesh slice
  type(sem_mesh_data) :: mesh_new
  integer :: iproc_new, nspec_new
  real(dp), allocatable :: model_gll_new(:,:,:,:,:)
  real(dp), allocatable :: xyz_center_new(:,:)

  !-- interpolation points
  real(dp), allocatable :: xyz_new(:,:)
  integer, allocatable :: idoubling_new(:)
  real(dp), parameter :: FILLVALUE_dp = huge(1.0_dp)
  integer :: igll, igllx, iglly, igllz, ngll_new
  real(dp), allocatable :: model_interp(:,:)

  !-- sem location
  type(sem_mesh_location), allocatable :: location_1slice(:)
  integer, parameter :: nnearest = 10
  real(dp) :: typical_size, max_search_dist, max_misloc, min_dist
  real(dp), allocatable :: misloc_final(:)
  integer, allocatable :: stat_final(:)

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
  read(args(1), '(a)') old_mesh_dir
  read(args(2), *) nproc_old
  read(args(3), '(a)') old_model_dir
  read(args(4), '(a)') model_tags
  read(args(5), '(a)') new_mesh_dir
  read(args(6), *) nproc_new
  read(args(7), '(a)') output_dir

  !===== parse model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)
  endif

  !===== interpolate old mesh/model onto new mesh

  ! typical element size at surface for old mesh
  ! read from values_from_mesher.h for old mesh
  typical_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                     ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                 * DEGREES_TO_RADIANS * R_UNIT_SPHERE

  max_search_dist = 10.0 * typical_size

  max_misloc = typical_size / 4.0

  !-- loop each new mesh slice

  do iproc_new = myrank, (nproc_new - 1), nrank

    print *, "# iproc_new=", iproc_new

    !-- read new mesh slice
    call sem_mesh_read(new_mesh_dir, iproc_new, iregion, mesh_new)

    !-- make arrays of xyz points
    !TODO use global points with different idoubling values instead of 
    ! all GLL points, as many GLL points are shared on the element faces
    nspec_new = mesh_new%nspec
    ngll_new = NGLLX * NGLLY * NGLLZ * nspec_new

    if (allocated(xyz_new)) then
      deallocate(xyz_new, idoubling_new, xyz_center_new)
    endif
    allocate(xyz_new(3, ngll_new), idoubling_new(ngll_new))
    allocate(xyz_center_new(3, nspec_new))

    do ispec = 1, nspec_new

      iglob = mesh_new%ibool(MIDX, MIDY, MIDZ, ispec)
      xyz_center_new(:, ispec) = mesh_new%xyz_glob(:, iglob)

      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            igll = igllx + &
                   NGLLX * ( (iglly-1) + &
                   NGLLY * ( (igllz-1) + &
                   NGLLZ * ( (ispec-1))))
            iglob = mesh_new%ibool(igllx, iglly, igllz, ispec)
            xyz_new(:, igll) = mesh_new%xyz_glob(:, iglob)
            idoubling_new(igll) = mesh_new%idoubling(ispec)
          enddo
        enddo
      enddo

    enddo

    !-- initialize variables for interpolation
    if (allocated(location_1slice)) then
      deallocate(location_1slice)
    endif
    allocate(location_1slice(ngll_new))

    if (allocated(stat_final)) then
      deallocate(stat_final, misloc_final)
    endif
    allocate(stat_final(ngll_new), misloc_final(ngll_new))
    stat_final = -1
    misloc_final = huge(1.0_dp)

    if (allocated(model_interp)) then
      deallocate(model_interp, model_gll_new)
    endif
    allocate(model_interp(nmodel, ngll_new))
    allocate(model_gll_new(nmodel, NGLLX, NGLLY, NGLLZ, nspec_new))
    model_interp = FILLVALUE_dp

    !-- read in the new model as the background model
!   call sem_io_read_gll_file_n(new_model_dir, iproc_new, iregion, &
!     model_names, nmodel, model_gll_new)
    model_gll_new = FILLVALUE_dp

    !-- loop each slices of the old mesh

    do iproc_old = 0, (nproc_old - 1)

      print *, "# iproc_old=", iproc_old

      ! read old mesh slice
      call sem_mesh_read(old_mesh_dir, iproc_old, iregion, mesh_old)

      nspec_old = mesh_old%nspec

      ! test if the new and old mesh slices are separated apart
      min_dist = huge(0.0_dp)
      do ispec = 1, nspec_old
        iglob = mesh_old%ibool(MIDX, MIDY, MIDZ, ispec)
        min_dist = min(min_dist, sqrt( minval( &
          (xyz_center_new(1,:) - mesh_old%xyz_glob(1,iglob))**2 + &
          (xyz_center_new(2,:) - mesh_old%xyz_glob(2,iglob))**2 + &
          (xyz_center_new(3,:) - mesh_old%xyz_glob(3,iglob))**2)))
      enddo
      if (min_dist > max_search_dist) then
        cycle
      endif

      ! read old model
      if (allocated(model_gll_old)) then
        deallocate(model_gll_old)
      endif
      allocate(model_gll_old(nmodel, NGLLX, NGLLY, NGLLZ, nspec_old))

      call sem_io_read_gll_file_n(old_model_dir, iproc_old, iregion, &
        model_names, nmodel, model_gll_old)

      ! locate points in this mesh slice
      call sem_mesh_locate_kdtree2(mesh_old, ngll_new, xyz_new, idoubling_new, &
        nnearest, max_search_dist, max_misloc, location_1slice)

      ! interpolate model only on points located inside an element
      do igll = 1, ngll_new

        ! skip points located inside more than one element 
        ! this would occur if the point lies exactly on the faces between two
        ! elements
        if (stat_final(igll)==1 .and. location_1slice(igll)%stat==1) then
!         print *, "[WARN] iproc_new=", iproc_new, &
!                  " iproc_old=", iproc_old, &
!                  " igll=", igll, &
!                  " this point is located inside more than one element!", &
!                  " some problem may occur.", &
!                  " only use the first located element."
          print *, "# multi-located ", xyz_new(:,igll)
          cycle
        endif

        ! for point located inside one element in the first time
        ! or closer to one element than located before
        if ( location_1slice(igll)%stat == 1 .or. &
            (location_1slice(igll)%stat == 0 .and. &
             location_1slice(igll)%misloc < misloc_final(igll)) ) &
        then

          ! interpolate model
          do imodel = 1, nmodel
            model_interp(imodel,igll) = &
              sum(location_1slice(igll)%lagrange * &
                model_gll_old(imodel, :, :, :, location_1slice(igll)%eid))
          enddo

          stat_final(igll)   = location_1slice(igll)%stat
          misloc_final(igll) = location_1slice(igll)%misloc

        endif

      enddo ! igll = 1, ngll_new

    enddo ! iproc_old

    !-- write out gll files for this new mesh slice

    ! reshape model_interp to model_gll
    do ispec = 1, nspec_new
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX

            igll = igllx + &
                   NGLLX * ( (iglly-1) + & 
                   NGLLY * ( (igllz-1) + & 
                   NGLLZ * ( (ispec-1))))

            if (stat_final(igll) /= -1) then
              model_gll_new(:, igllx, iglly, igllz, ispec) = &
                model_interp(:, igll)
            endif

          enddo
        enddo
      enddo
    enddo

    call sem_io_write_gll_file_n(output_dir, iproc_new, iregion, &
      model_names, nmodel, model_gll_new)

  enddo ! iproc_new

  !===== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
