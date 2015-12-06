subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_depth_mask"
  print '(a)', "    - make mask gll file (down-weighting shallow depth)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_source_depth_mask \"
  print '(a)', "    <nproc> <mesh_dir> <depth_stop> <depth_pass> <taper_type> "
  print '(a)', "    <out_tag> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (float) depth_stop: depth where the weight decreases to zero"
  print '(a)', "  (float) depth_pass: depth where the weight increases to one"
  print '(a)', "  (string) taper_type: taper function used between stop/pass deths (e.g. linear or cosine)"
  print '(a)', "  (string) out_tag: proc*_reg1_<out_tag>.bin"
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
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  real(dp) :: depth_stop
  real(dp) :: depth_pass
  character(len=MAX_STRING_LEN) :: taper_type
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
  ! depth mask
  integer :: taper_id
  real(dp) :: depth, depth_width

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
  read(args(3), *) depth_stop
  read(args(4), *) depth_pass 
  read(args(5), '(a)') taper_type
  read(args(6), '(a)') out_tag
  read(args(7), '(a)') out_dir

  ! validate inputs
  if (.not. (depth_pass>0.0 .and. depth_stop<depth_pass)) then
    if (myrank==0) then
      print *, "[ERROR] depth_pass > 0.0 and depth_stop < depth_pass"
    endif
    call abort_mpi()
  endif
  call synchronize_all()

  select case (trim(taper_type))
    case ('linear')
      taper_id = 1
    case ('cosine')
      taper_id = 2
    case default
      if (myrank==0) then
        print *, '[ERROR]: invalid taper_type'
        call abort_mpi()
      endif
  end select
  call synchronize_all()

  ! log output
  if (myrank == 0) then
    print *, "#[LOG] xsem_make_depth_mask"
    print *, "#[LOG] nproc=", nproc
    print *, "#[LOG] mesh_dir=", trim(mesh_dir)
    print *, "#[LOG] depth_stop=", depth_stop
    print *, "#[LOG] depth_pass=", depth_pass
    print *, "#[LOG] taper_type=", trim(taper_type)
    print *, "#[LOG] out_tag=", out_tag
    print *, "#[LOG] out_dir=", out_dir
  endif

  !===== loop each mesh slice

  ! non-dimensionalization
  depth_stop = depth_stop / R_EARTH_KM
  depth_pass = depth_pass / R_EARTH_KM
  depth_width = depth_pass - depth_stop

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

            ! depth mask
            if (depth < depth_stop) then
              weight = 0.0_dp
            elseif (depth > depth_pass) then
              weight = 1.0_dp
            else
              if (taper_id == 1) then
                weight = (depth - depth_stop)/depth_width
              elseif (taper_id == 2) then
                weight = 0.5 * (1.0 - cos(PI*(depth - depth_stop)/depth_width))
              else
                print *, "[ERROR] invalid taper_id"
                call abort_mpi()
              endif
            endif
           
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
