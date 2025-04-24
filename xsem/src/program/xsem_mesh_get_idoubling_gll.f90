subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_mesh_get_idoubling_gll"
  print '(a)', "    - make GLL files for element types (idoubling)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_mesh_get_idoubling_gll \"
  print '(a)', "    <nproc> <mesh_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_idoubling.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_mesh_get_idoubling_gll

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 3
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: idoubling_gll

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] check your inputs."
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') out_dir

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  allocate(idoubling_gll(NGLLX,NGLLY,NGLLZ,nspec))

  !====== loop all mesh slices
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! idoubling gll
    do ispec = 1, nspec
      idoubling_gll(:,:,:,ispec) = mesh_data%idoubling(ispec)
    enddo

    ! write out gll
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'idoubling', idoubling_gll)

  enddo

  !====== Finalize
  if (myrank == 0) close(IOUT)

  call synchronize_all()
  call finalize_mpi()

end program
