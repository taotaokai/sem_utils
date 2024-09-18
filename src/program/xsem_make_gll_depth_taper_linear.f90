subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_gll_depth_taper_linear "
  print '(a)', "    - make gll file of linear taper in depth"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_gll_depth_taper_linear \"
  print '(a)', "    <nproc> <mesh_dir> \"
  print '(a)', "    <depth1_km> <value1> <depth2_km> <value2> \"
  print '(a)', "    <out_tag> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (float) depth1_km: first depth from which gll value gradually changes from value1 to value2"
  print '(a)', "  (float) value1: gll value will be value1 from surface to depth1_km"
  print '(a)', "  (float) depth2_km: second depth at which the gll value reaches value2"
  print '(a)', "  (float) value2: gll value will be value2 below depth2_km"
  print '(a)', "  (string) out_tag: proc*_reg1_<out_tag>.bin"
  print '(a)', "  (string) out_dir: out_dir/proc*_reg1_<out_tag>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_make_kernel_gll

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
  real(dp) :: depth1, value1
  real(dp) :: depth2, value2
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
  ! depth taper gll
  real(dp), allocatable :: gll(:,:,:,:)
  real(dp) :: xyz(3), gllvalue, dvalue
  real(dp) :: depth, depth_width

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_make_kernel_gll: check your inputs."
      call abort_mpi()
    endif
  endif
  call synchronize_all()

  do i = 1, nargs
      call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), *) depth1
  read(args(4), *) value1 
  read(args(5), *) depth2
  read(args(6), *) value2
  read(args(7), '(a)') out_tag
  read(args(8), '(a)') out_dir

  ! validate inputs
  if (.not. (depth1>0.0 .and. depth2>depth1)) then
    if (myrank==0) then
      print *, "[ERROR] must satisfy: depth1 > 0.0 and depth1 < depth2"
    endif
    call abort_mpi()
  endif
  call synchronize_all()

  ! log output
  if (myrank == 0) then
    print *, "#[LOG] xsem_make_depth_gll"
    print *, "#[LOG] nproc=", nproc
    print *, "#[LOG] mesh_dir=", trim(mesh_dir)
    print *, "#[LOG] depth1,value1=", depth1, value1
    print *, "#[LOG] depth2,value2=", depth2, value2
    print *, "#[LOG] out_tag=", trim(out_tag)
    print *, "#[LOG] out_dir=", trim(out_dir)
  endif

  !===== loop each mesh slice

  ! non-dimensionalization
  depth1 = depth1 / EARTH_R_KM
  depth2 = depth2 / EARTH_R_KM
  depth_width = depth2 - depth1

  dvalue = value2 - value1

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    nspec = mesh_data%nspec

    ! initialize gll gll array 
    if (.not. allocated(gll)) then
      allocate(gll(NGLLX,NGLLY,NGLLZ,nspec))
    endif

    ! loop each gll point to set gll
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX

            ! gll point xyz
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)
            depth = 1.0 - sqrt(sum(xyz**2))

            ! depth gll
            if (depth < depth1) then
              gllvalue = value1
            elseif (depth > depth2) then
              gllvalue = value2
            else
              gllvalue = value1 + dvalue*(depth - depth1)/depth_width
            endif
           
            gll(igllx,iglly,igllz,ispec) = gllvalue

          enddo
        enddo
      enddo
    enddo

    print *,'min/max=', minval(gll), maxval(gll)

    ! save gll gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_tag, gll)

  enddo ! iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
