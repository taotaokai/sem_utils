subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_source_depth_mask "
  print '(a)', "    - make mask gll file (down-weighting source region and shallow depth)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_source_mask "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <source_xyz_list> <source_gaussa> "
  print '(a)', "    <surface_weight> <depth_gaussa> "
  print '(a)', "    <out_dir> <out_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) source_xyz_list: list of source locations(x,y,z) in SEM (normalized)"
  print '(a)', "  (float) source_gaussa: Gaussian width (one sigma) in km"
  print '(a)', "  (float) surface_weight: weighting at the surface, should be between (0, 1), usually small number, say 0.001"
  print '(a)', "  (float) depth_gaussa: Gaussian width (one sigma) in km"
  print '(a)', "  (string) out_dir: <out_dir>/proc*_reg1_<out_name>.bin"
  print '(a)', "  (string) out_name: file name will be proc*_reg1_<out_name>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_make_source_depth_mask

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: source_xyz_list
  real(dp) :: source_gaussa
  real(dp) :: surface_weight
  real(dp) :: depth_gaussa
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank
  ! source locations
  character(len=MAX_STRING_LEN), allocatable :: lines(:)
  integer :: nsource, isrc
  real(dp), allocatable :: source_xyz(:,:)
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec, iglob, igllx, iglly, igllz, ispec
  ! mask gll 
  real(dp), allocatable :: mask(:,:,:,:)
  real(dp) :: xyz(3), weight
  ! source weighting
  real(dp) :: dist_sq
  ! depth weighting
  real(dp) :: depth

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_make_source_mask: check your inputs."
      call abort_mpi()
    endif
  endif
  call synchronize_all()

  do i = 1, nargs
      call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') source_xyz_list
  read(args(4), *) source_gaussa
  read(args(5), *) surface_weight
  read(args(6), *) depth_gaussa
  read(args(7), '(a)') out_dir
  read(args(8), '(a)') out_name

  ! validate input
  if (myrank == 0) then
    if (source_gaussa < 0.0 .or. depth_gaussa < 0.0) then
      call selfdoc()
      print *, "[ERROR] gaussa must > 0.0"
      call abort_mpi()
    endif
    if (surface_weight <= 0.0 .or. surface_weight >= 1.0) then
      call selfdoc()
      print *, "[ERROR] surface_weight must between (0,1)"
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  ! log output
  if (myrank == 0) then
    print *, "#[LOG] xsem_make_source_depth_mask"
    print *, "#[LOG] nproc=", nproc
    print *, "#[LOG] mesh_dir=", trim(mesh_dir)
    print *, "#[LOG] source_xyz_list=", trim(source_xyz_list)
    print *, "#[LOG] source_gaussa=", source_gaussa
    print *, "#[LOG] surface_weight=", surface_weight
    print *, "#[LOG] depth_gaussa=", depth_gaussa
    print *, "#[LOG] out_name=", trim(out_name)
    print *, "#[LOG] out_dir=", trim(out_dir)
  endif

  !====== read source_xyz_list
  call sem_utils_read_line(source_xyz_list, lines, nsource)

  allocate(source_xyz(3, nsource))
  do isrc = 1, nsource
    ! read source x,y,z
    read(lines(isrc), *) source_xyz(1, isrc), source_xyz(2, isrc), source_xyz(3, isrc)
  enddo

  !===== loop each mesh slice

  ! non-dimensionalization
  source_gaussa = source_gaussa / R_EARTH_KM
  depth_gaussa = depth_gaussa / R_EARTH_KM

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    nspec = mesh_data%nspec

    ! initialize mask gll array 
    if (.not. allocated(mask)) then
      allocate(mask(NGLLX,NGLLY,NGLLZ,nspec))
    endif

    ! loop each gll point to set mask
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX

            ! gll point xyz
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            ! source mask: mask source region
            weight = 1.0
            do isrc = 1, nsource
              dist_sq = sum((xyz - source_xyz(:,isrc))**2)
              weight = weight * (1.0 - exp(-0.5*dist_sq/source_gaussa**2))
            enddo

            ! depth weighting
            depth = 1.0 - sqrt(sum(xyz**2))
            if (depth < 0.0 ) then
              depth = 0.0
            endif
            weight = weight * (surface_weight + (1.0 - surface_weight)*(1.0 - exp(-0.5*(depth/depth_gaussa)**2)))

            ! mask gll
            mask(igllx,iglly,igllz,ispec) = weight

          enddo
        enddo
      enddo
    enddo

    print *,'min/max=', minval(mask), maxval(mask)

    ! save mask gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, mask)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
