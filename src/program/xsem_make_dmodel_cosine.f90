subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_make_dmodel_cosine"
  print '(a)', "    - make plane wave like model perturbation"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_make_dmodel_random "
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <amplitude> <wave_length_km> <nx> <ny> <nz>"
  print '(a)', "    <out_dir> <out_name> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', "  dmodel = A0 * cos(2*PI/wave_length*<n, x>)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (float) amplitude:  A0 "
  print '(a)', "  (float) wave_length_km:  wave length in kilometers"
  print '(a)', "  (float) nx:  wave vector along x axis "
  print '(a)', "  (float) ny:  wave vector along y axis "
  print '(a)', "  (float) nz:  wave vector along z axis "
  print '(a)', "  (string) out_dir:  out directory for dmodel"
  print '(a)', "  (string) out_name:  proc000***_reg1_<out_name>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_make_dmodel_cosine

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 9
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  real(dp) :: amplitude
  real(dp) :: wave_length
  real(dp) :: nx
  real(dp) :: ny
  real(dp) :: nz 
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec, ispec, igllz, iglly, igllx, iglob
  real(dp) :: xyz(3), norm, wave_number
  ! dmodel
  real(dp), dimension(:,:,:,:), allocatable :: dmodel

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] wrong number of inputs!"
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),*) amplitude 
  read(args(4),*) wave_length
  read(args(5),*) nx 
  read(args(6),*) ny 
  read(args(7),*) nz 
  read(args(8),'(a)') out_dir
  read(args(9),'(a)') out_name

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(dmodel(NGLLX,NGLLY,NGLLZ,nspec))

  ! non-dimensionalize 
  wave_number = TWO_PI/(wave_length/R_EARTH_KM)
  norm = (nx**2 + ny**2 + nz**2)**0.5
  nx = nx/norm
  ny = ny/norm
  nz = nz/norm
 
  !====== create new model
  do iproc = myrank, (nproc-1), nrank

    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! create cosine plan wave
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            dmodel(igllx,iglly,igllz,ispec) = amplitude*cos(wave_number * (nx*xyz(1)+ny*xyz(2)+nz*xyz(3)))

          enddo
        enddo
      enddo
    enddo

    ! write out dmodel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, dmodel)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
