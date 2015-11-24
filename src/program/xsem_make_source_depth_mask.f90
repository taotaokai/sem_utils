subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_source_depth_mask "
  print '(a)', "    - make mask gll file (mask source, shallow depth)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_kernel_mask \"
  print '(a)', "    <nproc> <mesh_dir> <source_xyz_list> <source_mask_radius> \"
  print '(a)', "    <stop_depth> <pass_depth> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) source_xyz_list: list of source locations(x,y,z) in SEM (normalized)"
  print '(a)', "  (float) source_mask_radius: Gaussian width (one sigma) in km"
  print '(a)', "  (float) stop_depth: stop depth (km) where mask = 0"
  print '(a)', "  (float) pass_depth: pass depth (km) where mask = 1"
  print '(a)', "  (string) out_dir: output directory for proc*_reg1_source_mask.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_make_kernel_mask

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: source_xyz_list
  real(dp) :: source_mask_radius
  real(dp) :: stop_depth 
  real(dp) :: pass_depth 
  character(len=MAX_STRING_LEN) :: out_dir

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
  ! source
  real(dp) :: dist_sq
  ! depth mask
  real(dp) :: depth

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_make_kernel_mask: check your inputs."
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
  read(args(4), *) source_mask_radius
  read(args(5), *) stop_depth
  read(args(6), *) pass_depth
  read(args(7), '(a)') out_dir

  !====== read source_xyz_list
  call sem_utils_read_line(source_xyz_list, lines, nsource)

  allocate(source_xyz(3, nsource))
  do isrc = 1, nsource
    ! read source x,y,z 
    read(lines(isrc), *) source_xyz(1, isrc), source_xyz(2, isrc), &
                         source_xyz(3, isrc)
  enddo

  !===== loop each mesh slice

  !non-dimensionalize
  source_mask_radius = source_mask_radius / R_EARTH_KM
  stop_depth = stop_depth / R_EARTH_KM
  pass_depth = pass_depth / R_EARTH_KM

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
            depth = 1.0 - sqrt(sum(xyz**2))
            
            weight = 1.0_dp
            ! source mask: mask source region
            if (source_mask_radius > 0.0) then
              do isrc = 1, nsource
                dist_sq = sum((xyz - source_xyz(:,isrc))**2)
                weight = weight * (1.0_dp - &
                    exp(-0.5*dist_sq/source_mask_radius**2))
              enddo
            endif
            ! depth mask: mask shallow depth
            if (pass_depth > stop_depth) then
              if (stop_depth<depth .and. depth<pass_depth) then
                weight = weight * (0.5 - 0.5*cos(PI* & 
                  (depth-stop_depth)/(pass_depth-stop_depth)))
              elseif (depth <= stop_depth) then
                weight = 0.0_dp
              endif
            endif

            mask(igllx,iglly,igllz,ispec) = weight

          enddo
        enddo
      enddo
    enddo

    print *,'mask value range: ', minval(mask), maxval(mask)

    ! save mask gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, "mask", mask)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
