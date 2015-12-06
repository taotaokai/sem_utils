subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_source_mask "
  print '(a)', "    - make mask gll file (down-weighting around source region)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_source_mask \"
  print '(a)', "    <nproc> <mesh_dir> <source_xyz_list> <source_gaussa> \"
  print '(a)', "    <out_tag> <out_dir> "
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
  print '(a)', "  (string) out_tag: tag in file name proc*_reg1_<out_tag>.bin"
  print '(a)', "  (string) out_dir: out_dir/proc*_reg1_<out_tag>.bin"
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
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: source_xyz_list
  real(dp) :: source_gaussa
  character(len=MAX_STRING_LEN) :: out_tag
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec, iglob, igllx, iglly, igllz, ispec
  ! mask gll 
  real(dp), allocatable :: mask(:,:,:,:)
  real(dp) :: xyz(3), weight
  ! source
  real(dp) :: dist_sq

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
  read(args(5), '(a)') out_tag
  read(args(6), '(a)') out_dir

  ! validate input
  if (source_gaussa < 0.0) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] source_gaussa must > 0.0"
      call abort_mpi()
    endif
  endif
  call synchronize_all()

  ! log output
  if (myrank == 0) then
    print *, "#[LOG] xsem_make_source_mask"
    print *, "#[LOG] nproc=", nproc
    print *, "#[LOG] mesh_dir=", trim(mesh_dir)
    print *, "#[LOG] source_xyz_list=", trim(source_xyz_list)
    print *, "#[LOG] source_gaussa=", source_gaussa
    print *, "#[LOG] out_tag=", out_tag
    print *, "#[LOG] out_dir=", out_dir
  endif

  !====== read source_xyz_list
  call sem_utils_read_line(source_xyz_list, lines, nsource)

  allocate(source_xyz(3, nsource))
  do isrc = 1, nsource
    ! read source x,y,z 
    read(lines(isrc), *) source_xyz(1, isrc), source_xyz(2, isrc), &
                         source_xyz(3, isrc)
  enddo

  !===== loop each mesh slice

  ! non-dimensionalization
  source_gaussa = source_gaussa / R_EARTH_KM

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
            do isrc = 1, nsource
              dist_sq = sum((xyz - source_xyz(:,isrc))**2)
              weight = weight * (1.0_dp - &
                  exp(-0.5*dist_sq/source_gaussa**2))
            enddo

            mask(igllx,iglly,igllz,ispec) = weight

          enddo
        enddo
      enddo
    enddo

    print *,'min/max=', minval(mask), maxval(mask)

    ! save mask gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_tag, mask)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
