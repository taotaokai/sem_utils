subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_crust_to_1d_gll"
  print '(a)', "    add crustal elements (idoubling==IFLAG_CRUST) into a 1D gll model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_crust_to_1d_gll \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir_with_crust> <model_dir_1d> <model_name> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir_with_crust:  model files with crustal elements"
  print '(a)', "  (string) model_dir_1d:  model files of only 1d structure"
  print '(a)', "  (string) model_name: e.g. vsv"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_<model_name>"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_add_crust_to_1d_gll

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir_with_crust
  character(len=MAX_STRING_LEN) :: model_dir_1d
  character(len=MAX_STRING_LEN) :: model_name
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! IFLAG
  integer :: iflag, idoubling
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: model_new, model_with_crust, model_1d

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] check your inputs."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') model_dir_with_crust
  read(args(4),'(a)') model_dir_1d
  read(args(5),'(a)') model_name 
  read(args(6),'(a)') out_dir

  !====== get mesh geometry
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! intialize arrays 
  allocate(model_with_crust(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(model_1d(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate thomsen parameters
  do iproc = 0, (nproc-1)
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! read models
    call sem_io_read_gll_file_1(model_dir_with_crust, iproc, iregion, model_name, model_with_crust)
    call sem_io_read_gll_file_1(model_dir_1d, iproc, iregion, model_name, model_1d)

    do ispec = 1, nspec
      idoubling = mesh_data%idoubling(ispec)
      ! use old IFLAG in case of mesh layer re-number (as done in sem_mesh_read)
      if (   idoubling == (10*IFLAG_CRUST + 0) &
        .or. idoubling == (10*IFLAG_CRUST + 1) &
        .or. idoubling == (10*IFLAG_CRUST + 2) ) then
        idoubling = IFLAG_CRUST ! 3 layer crustal mesh
      endif
      if (idoubling == IFLAG_CRUST) then
        model_1d(:,:,:,ispec) = model_with_crust(:,:,:,ispec)
      endif
    enddo

    ! write new model
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, model_name, model_1d)

  enddo

end program
