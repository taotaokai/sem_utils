subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_interp_xyz - interpolate SEM model on given points "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_interp_xyz \ "
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <xyz_list> <out_list>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  interpolate SEM mesh and associated model values onto the given xyz list"
  print '(a)', "    the SEM GLL interpolation function is used."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices "
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (string) xyz_list:  input file name, "
  print '(a)', "    each row: <x> <y> <z>"
  print '(a)', "  (string) out_list:  output file name, "
  print '(a)', "    each row: <x> <y> <z> misloc <misloc> [<tag> <val>]... "
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_interp_xyz
  
  use sem_constants
  use sem_parallel
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  !-- command line args
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags
  integer :: nproc
  character(len=MAX_STRING_LEN) :: xyz_list, out_list

  !-- local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier 

  !-- model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- interpolation points
  integer :: ipoint, npoint
  character(len=MAX_STRING_LEN), allocatable :: xyz_strings(:) 
  real(dp), allocatable :: xyz(:,:)
  integer, allocatable :: idoubling(:)

  !-- sem location 
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec, iglob
  real(dp) :: dist
  type(sem_mesh_location), allocatable :: location_1slice(:)
  integer, parameter :: nnearest = 10
  real(dp) :: typical_size, max_search_dist, max_misloc
  real(dp), allocatable :: misloc_myrank(:), misloc_gather(:), misloc_copy(:)
  integer, allocatable :: stat_myrank(:), stat_gather(:), stat_copy(:)
  !-- model interpolation
  real(sp), parameter :: FILLVALUE_sp = huge(0.0_sp)
  real(dp), allocatable :: model_gll(:,:,:,:,:)
  real(dp), allocatable :: model_interp_myrank(:,:)
  real(dp), allocatable :: model_interp_gather(:,:)
  real(dp), allocatable :: model_interp_copy(:,:)

  !-- mpi 
  integer :: myrank, nrank, iworker
  ! mpi_send/recv
  integer, parameter :: MPI_TAG_stat = 10
  integer, parameter :: MPI_TAG_misloc = 11
  integer, parameter :: MPI_TAG_model_interp = 12

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_interp_xyz: check your input arguments."
      call abort_mpi()
      stop
    endif
  endif
   
  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), '(a)') xyz_list
  read(args(6), '(a)') out_list

  !===== read xyz list

  call sem_utils_read_line(xyz_list, xyz_strings, npoint)

  allocate(xyz(3, npoint))

  do i = 1, npoint
    read(xyz_strings(i),*,iostat=ier) xyz(1,i), xyz(2,i), xyz(3,i)
    if (ier /= 0) then
      print *, "[ERROR] failed to read (x y z) value at line ", i
      stop
    endif
  enddo

  print *, '# npoint=', npoint

  ! layer ID
  allocate(idoubling(npoint))
  idoubling = IFLAG_DUMMY

  call synchronize_all()

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)
  endif

  call synchronize_all()

  !===== locate xyz in each mesh slice

  !-- initialize variables
  allocate(location_1slice(npoint))

  allocate(stat_myrank(npoint), misloc_myrank(npoint))
  stat_myrank = -1
  misloc_myrank = HUGE(1.0_dp)

  allocate(model_interp_myrank(nmodel, npoint))
  model_interp_myrank = FILLVALUE_sp

  ! typical element size at surface
  typical_size = max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                     ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                 * DEGREES_TO_RADIANS * R_UNIT_SPHERE

  max_search_dist = 5.0 * typical_size

  max_misloc = typical_size / 4.0

  !-- loop mesh slices for myrank
  loop_proc: do iproc = myrank, (nproc - 1), nrank

    print *, "# iproc=", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! test if the interpolation points and mesh are separated apart
    ! only use mesh element centers
    dist = huge(0.0_dp)
    do ispec = 1, nspec
      iglob = mesh_data%ibool(MIDX, MIDY, MIDZ, ispec)
      dist = min(dist, sqrt( minval( &
        (xyz(1,:) - mesh_data%xyz_glob(1,iglob))**2 + &
        (xyz(2,:) - mesh_data%xyz_glob(2,iglob))**2 + &
        (xyz(3,:) - mesh_data%xyz_glob(3,iglob))**2)))
    enddo
    if (dist > max_search_dist) then
      cycle
    endif

    ! read model
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(nmodel, NGLLX, NGLLY, NGLLZ, nspec))

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
      model_names, nmodel, model_gll)

    ! locate points in this mesh slice
    call sem_mesh_locate_kdtree2(mesh_data, npoint, xyz, idoubling, &
      nnearest, max_search_dist, max_misloc, location_1slice)

    print *, "# iproc=", iproc, &
      " number of located points: ", count(location_1slice%stat>=0)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! safety check
      if (stat_myrank(ipoint) == 1 .and. &
          location_1slice(ipoint)%stat == 1 ) then
        print *, "[WARN] iproc=", iproc, " ipoint=", ipoint, &
          " This point is located inside more than one element!", &
          " Only use the first located element."
        cycle
      endif

      ! for point located inside one element but not happened before 
      ! or closer to one element than located before
      if (location_1slice(ipoint)%stat == 1 &
          .or. &
          (location_1slice(ipoint)%stat == 0 .and. &
           location_1slice(ipoint)%misloc < misloc_myrank(ipoint)) ) &
      then

        ! interpolate model
        do imodel = 1, nmodel
          model_interp_myrank(imodel, ipoint) = &
            sum( location_1slice(ipoint)%lagrange * &
                 model_gll(imodel, :, :, :, location_1slice(ipoint)%eid) )
        enddo

        stat_myrank(ipoint)   = location_1slice(ipoint)%stat
        misloc_myrank(ipoint) = location_1slice(ipoint)%misloc

      endif

    enddo ! ipoint

  enddo loop_proc

  !-- gather all location results from work processors
  call synchronize_all()

  if (myrank /= 0) then

    ! send stat/misloc/model_interp on myrank to MASTER 

    call send_i(stat_myrank, npoint, 0, MPI_TAG_stat)
    call send_dp(misloc_myrank, npoint, 0, MPI_TAG_misloc)
    call send2_dp(model_interp_myrank, nmodel*npoint, 0, MPI_TAG_model_interp)

  else ! gather data on MASTER processor (0)

    allocate(stat_gather(npoint), misloc_gather(npoint), &
      model_interp_gather(nmodel, npoint))

    allocate(stat_copy(npoint), misloc_copy(npoint), &
      model_interp_copy(nmodel, npoint))

    stat_gather = stat_myrank 
    misloc_gather = misloc_myrank 
    model_interp_gather = model_interp_myrank 

    do iworker = 1, (nrank -1)

      ! receive stat/misloc/model_interp from each worker processes 
      call recv_i(stat_copy, npoint, iworker , MPI_TAG_stat)
      call recv_dp(misloc_copy, npoint, iworker, MPI_TAG_misloc)
      call recv2_dp(model_interp_copy, nmodel*npoint, iworker, &
        MPI_TAG_model_interp)

      ! gather data
      do ipoint = 1, npoint

        ! safety check
        if (stat_gather(ipoint) == 1) then
          if (stat_copy(ipoint) == 1) then
            print *, "[WARN] ipoint=", ipoint
            print *, "------ this point is located inside more than one element!"
            print *, "------ some problem may occur."
            print *, "------ only use the first located element."
          endif
          cycle
        endif

        ! for point located inside one element but not happened before 
        ! or closer to one element than located before
        if (stat_copy(ipoint) == 1 &
            .or. &
            (stat_copy(ipoint) == 0 .and. &
             misloc_copy(ipoint) < misloc_gather(ipoint)) ) &
        then
          stat_gather(ipoint) = stat_copy(ipoint)
          misloc_gather(ipoint) = misloc_copy(ipoint)
          model_interp_gather(:, ipoint) = model_interp_copy(:,ipoint) 
        endif

      enddo ! ipoint

    enddo ! iworker

    ! convert misloc relative to typical element size
    where (stat_gather /= -1)
      misloc_gather = misloc_gather / typical_size
    endwhere

  endif ! myrank /= 0

  call synchronize_all()

  !===== write out results on the master process 

  if (myrank == 0) then

    open(unit=IOUT, file=out_list, status='unknown', iostat=ier)
    if (ier /= 0) then
      print *, '[ERROR] failed to open file: ', out_list
      call abort_mpi()
      stop
    endif

    do ipoint = 1, npoint

      write(IOUT, '(3E15.7,a,E15.7)', advance='no') &
        xyz(:,ipoint),' misloc ', misloc_gather(ipoint)

      do imodel = 1, nmodel
        write(IOUT, '(2x,a,2x,E15.7)', advance='no') &
          trim(model_names(imodel)), model_interp_gather(imodel, ipoint)
      enddo

      write(IOUT,'(a)') ""

    enddo

    close(IOUT)

  endif

  !===== Finalize MPI

  call synchronize_all()
  call finalize_mpi()

end program
