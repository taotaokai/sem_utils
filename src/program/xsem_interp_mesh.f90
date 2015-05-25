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
! print '(a)', "    <max_search_dist> <max_misloc> "
  print '(a)', "    <new_mesh_dir> <nproc_new> <new_model_dir> "
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
! print '(a)', "  (float) max_search_dist: maximum distance to "
! print '(a)', "    the n nearest elements from target point "
! print '(a)', "  (float) max_misloc: maximum location error for points "
! print '(a)', "    lying outside target mesh/layer volume "
  print '(a)', "  (string) new_mesh_dir: directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (int) nproc_new: number of slices of the new mesh"
  print '(a)', "  (string) new_model_dir: output directory for new model files"
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  1) This program must run in parallel, e.g. mpirun -n <nproc> ..."
  print '(a)', "  2) use values_from_mesh.h for old mesh to compile"
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

  character(len=MAX_STRING_LEN) :: new_mesh_dir, new_model_dir
  integer :: nproc_new

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
  integer :: iproc_new, nspec_new, nglob_new
  real(dp), allocatable :: model_gll_new(:,:,:,:,:)

  !-- interpolation points
  real(dp), allocatable :: xyz_interp(:,:)
  integer, allocatable :: idoubling_interp(:)
  integer, allocatable :: ispec_interp(:), iglob_interp(:)
  integer, allocatable :: RANGE_1_NPOINT(:)
  integer, parameter :: nshare_max = 2
  integer :: idoubling_share(nshare_max)
  integer :: ncount
  integer :: nshare, ipoint, npoint_interp, idoubling1
  real(dp), parameter :: FILLVALUE_dp = huge(1.0_dp)
  real(dp), allocatable :: model_interp(:,:)

  !-- sem location
  type(sem_mesh_location), allocatable :: location_1slice(:)
  integer, parameter :: nnearest = 10
  real(dp) :: typical_size, max_search_dist, max_misloc
  real(dp), allocatable :: misloc_final(:)
  integer, allocatable :: stat_final(:)

  !-- model output
  logical, allocatable :: idx_ispec(:), idx(:)
  integer :: igllx, iglly, igllz
  integer :: ipoint1(1)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      stop "[ERROR] xsem_interp_mesh: check your input arguments."
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
  read(args(7), '(a)') new_model_dir

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

  max_search_dist = 5.0 * typical_size

  max_misloc = typical_size / 4.0

  !-- loop each new mesh slice

  do iproc_new = myrank, (nproc_new - 1), nrank

    print *, "# iproc_new=", iproc_new

    !-- read new mesh slice
    call sem_mesh_read(new_mesh_dir, iproc_new, iregion, mesh_new)

    nspec_new = mesh_new%nspec
    nglob_new = mesh_new%nglob

    !-- make arrays of interpolation points

    ! count the number of interpolation points from global points
    ! the same global point with different idoubling values is taken 
    ! as different interpolation points
    npoint_interp = 0
    do iglob = 1, nglob_new

      ! count number of different idoubling values shared at one global point
      nshare = 0
      idoubling_share = huge(0) !make sure this value is NOT used as any layer id's
      do ispec = 1, nspec_new
        ncount = count(mesh_new%ibool(:,:,:,ispec) == iglob)
        ! one global point can only occur once in each element
        if (ncount == 1) then
          idoubling1 = mesh_new%idoubling(ispec)
          ! store idoubling1 if it is unique
          if (all(idoubling_share /= idoubling1)) then
            nshare = nshare + 1
            if (nshare > nshare_max) then
              print *, "[ERROR] iglob=", iglob, " ispec=", ispec, &
                       " nshare is greater than nshare_max"
              call abort_mpi()
            endif
            idoubling_share(nshare) = idoubling1
          endif
        elseif (ncount /= 0) then
          print *, "[ERROR] ncount=", ncount, " ispec=", ispec, &
                   " iglob=", iglob, &
                   " two or more GLL points in one element refer to the same", &
                   " global point. The mesh file could be wrong."
          call abort_mpi()
        endif
      enddo

      npoint_interp = npoint_interp + nshare

    enddo

    ! initialize arrays 
    if (allocated(xyz_interp)) then
      deallocate(xyz_interp, idoubling_interp, ispec_interp)
    endif
    allocate(xyz_interp(3, npoint_interp), idoubling_interp(npoint_interp))
    allocate(ispec_interp(npoint_interp), iglob_interp(npoint_interp))
    allocate(RANGE_1_NPOINT(npoint_interp))
    RANGE_1_NPOINT = [ (ipoint, ipoint=1,npoint_interp) ]

    ! now get interpolation points
    ! TODO: anyway to avoid the dumplication?
    ipoint = 0
    do iglob = 1, nglob_new
      ! count number of different idoubling values shared at one global point
      nshare = 0
      idoubling_share = huge(0) ! make sure this value is NOT used as any layer id's
      do ispec = 1, nspec_new
        ncount = count(mesh_new%ibool(:,:,:,ispec) == iglob)
        ! one global point can only occur once in each element
        if (ncount == 1) then
          idoubling1 = mesh_new%idoubling(ispec)
          ! store idoubling1 if it is unique
          if (all(idoubling_share /= idoubling1)) then
            nshare = nshare + 1
            idoubling_share(nshare) = idoubling1
            ! interpolation points
            ipoint = ipoint + 1
            xyz_interp(:, ipoint) = mesh_new%xyz_glob(:, iglob)
            idoubling_interp(ipoint) = mesh_new%idoubling(ispec)
            iglob_interp(ipoint) = iglob
            ispec_interp(ipoint) = ispec
          endif
        endif
      enddo
    enddo

!   npoint_interp = NGLLX * NGLLY * NGLLZ * nspec_new
!   do ispec = 1, nspec_new
!     do igllz = 1, NGLLZ
!       do iglly = 1, NGLLY
!         do igllx = 1, NGLLX
!           igll = igllx + &
!                  NGLLX * ( (iglly-1) + &
!                  NGLLY * ( (igllz-1) + &
!                  NGLLZ * ( (ispec-1))))
!           iglob = mesh_new%ibool(igllx, iglly, igllz, ispec)
!           xyz_interp(:, igll) = mesh_new%xyz_glob(:, iglob)
!           idoubling_interp(igll) = mesh_new%idoubling(ispec)
!         enddo
!       enddo
!     enddo
!   enddo

    !-- initialize variables for interpolation
    if (allocated(location_1slice)) then
      deallocate(location_1slice)
    endif
    allocate(location_1slice(npoint_interp))

    if (allocated(stat_final)) then
      deallocate(stat_final, misloc_final)
    endif
    allocate(stat_final(npoint_interp), misloc_final(npoint_interp))
    stat_final = -1
    misloc_final = huge(1.0_dp)

    if (allocated(model_interp)) then
      deallocate(model_interp, model_gll_new)
    endif
    allocate(model_interp(nmodel, npoint_interp))
    allocate(model_gll_new(nmodel, NGLLX, NGLLY, NGLLZ, nspec_new))
    model_interp = FILLVALUE_dp
 
    !-- loop each slices of the old mesh

    do iproc_old = 0, (nproc_old - 1)

      print *, "# iproc_old=", iproc_old

      ! read old mesh slice
      call sem_mesh_read(old_mesh_dir, iproc_old, iregion, mesh_old)

      ! read old model
      nspec_old = mesh_old%nspec
      if (allocated(model_gll_old)) then
        deallocate(model_gll_old)
      endif
      allocate(model_gll_old(nmodel, NGLLX, NGLLY, NGLLZ, nspec_old))

      call sem_io_read_gll_file_n(old_model_dir, iproc_old, iregion, &
        model_names, nmodel, model_gll_old)

      ! locate points in this mesh slice
      call sem_mesh_locate_kdtree2(mesh_old, npoint_interp, xyz_interp, idoubling_interp, &
        nnearest, max_search_dist, max_misloc, location_1slice)

      ! interpolate model only on points located inside an element
      do ipoint = 1, npoint_interp

        ! skip points located inside more than one element 
        ! this would occur if the point lies exactly on the faces between two
        ! elements
        if (stat_final(ipoint)==1 .and. location_1slice(ipoint)%stat==1) then
!         print *, "[WARN] iproc_new=", iproc_new, &
!                  " iproc_old=", iproc_old, &
!                  " ipoint=", ipoint, &
!                  " this point is located inside more than one element!", &
!                  " some problem may occur.", &
!                  " only use the first located element."
          print *, "# multi-located ", xyz_interp(:,ipoint)
          cycle
        endif

        ! for point located inside one element in the first time
        ! or closer to one element than located before
        if ( location_1slice(ipoint)%stat == 1 .or. &
             (location_1slice(ipoint)%stat == 0 .and. &
              location_1slice(ipoint)%misloc < misloc_final(ipoint)) ) &
        then

          ! interpolate model
          do imodel = 1, nmodel
            model_interp(imodel,ipoint) = &
              sum(location_1slice(ipoint)%lagrange * &
                model_gll_old(imodel, :, :, :, location_1slice(ipoint)%eid))
          enddo

          stat_final(ipoint)   = location_1slice(ipoint)%stat
          misloc_final(ipoint) = location_1slice(ipoint)%misloc

        endif

      enddo ! ipoint = 1, npoint_interp

    enddo ! iproc_old

    !-- write out gll files for this new mesh slice

    if (allocated(idx_ispec)) then
      deallocate(idx_ispec)
    endif
    allocate(idx_ispec(npoint_interp), idx(npoint_interp))

    ! put model_interp back into model_gll shape
    do ispec = 1, nspec_new

      idx_ispec = (ispec_interp == ispec)

      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX

            iglob = mesh_new%ibool(igllx, iglly, igllz, ispec)

            idx = idx_ispec .and. (iglob_interp == iglob)

            ncount = count(idx)
            if (ncount /= 1) then
              print *, "[ERROR] ispec_new=", ispec, " ncount=", ncount, &
                " two or more interpolation points foound for one GLL point"
              call abort_mpi()
              stop
            endif

            ipoint1 = pack(RANGE_1_NPOINT, mask=idx)
            model_gll_new(:, igllx, iglly, igllz, ispec) = &
              model_interp(:, ipoint1(1))

          enddo
        enddo
      enddo
    enddo

!   ! reshape model_interp to model_gll
!   do ispec = 1, nspec_new
!     do igllz = 1, NGLLZ
!       do iglly = 1, NGLLY
!         do igllx = 1, NGLLX
!           igll = igllx + &
!                  NGLLX * ( (iglly-1) + & 
!                  NGLLY * ( (igllz-1) + & 
!                  NGLLZ * ( (ispec-1))))
!           if (stat_final(igll) /= -1) then
!             model_gll_new(:, igllx, iglly, igllz, ispec) = &
!               model_interp(:, igll)
!           endif
!         enddo
!       enddo
!     enddo
!   enddo

    call sem_io_write_gll_file_n(new_model_dir, iproc_new, iregion, &
      model_names, nmodel, model_gll_new)

  enddo ! iproc_new

  !===== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
