subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_tiso_enforce_iso_iflag"
  print '(a)', "    enforce isotropic element in a vti model (vpv,vph,vsv,vsh,eta) based on model iflag"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_model_tiso_enforce_iso_iflag \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <flag_str> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vpv,vph,vsv,vsh,eta].bin"
  print '(a)', "  (string) flag_str:  element flag e.g. crust,moho_80,80_220,220_670,lower_mantle"
  print '(a)', "  (string) out_dir:  output directory for proc*_reg1_[eps,gamma,delta,eps_minus_delta]"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
end subroutine


program xsem_model_tiso_enforce_iso_iflag

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: flag_str
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
  real(dp), dimension(:,:,:,:), allocatable :: vph, vpv, vsv, vsh, eta

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
  read(args(3),'(a)') model_dir 
  read(args(4),'(a)') flag_str
  read(args(5),'(a)') out_dir

  call synchronize_all()

  !====== get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)
  call synchronize_all()

  ! intialize arrays 
  allocate(vph(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta(NGLLX,NGLLY,NGLLZ,nspec))

  !====== parse input flags
  iflag = -1
  if (trim(flag_str) == 'crust') then
    iflag = IFLAG_CRUST
  else if (trim(flag_str) == 'moho_80') then
    iflag = IFLAG_80_MOHO
  else if (trim(flag_str) == '80_220') then
    iflag = IFLAG_220_80
  else if (trim(flag_str) == '220_670') then
    iflag = IFLAG_670_220
  else if (trim(flag_str) == 'lower_mantle') then
    iflag = IFLAG_MANTLE_NORMAL
  else
    print *, "[ERROR] unrecognized IFLAG type: ", trim(flag_str)
    stop
  endif

  !====== calculate thomsen parameters
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read mesh
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)

    ! enforce isotropic crust with IFLAG_CRUST, 
    ! FIXME incorrect with moho stretching
    do ispec = 1, nspec
      idoubling = mesh_data%idoubling(ispec)
      ! use old IFLAG in case of mesh layer re-number (as done in sem_mesh_read)
      if (   idoubling == (10*IFLAG_CRUST + 0) &
        .or. idoubling == (10*IFLAG_CRUST + 1) &
        .or. idoubling == (10*IFLAG_CRUST + 2) ) then
        idoubling = IFLAG_CRUST ! 3 layer crustal mesh
      endif
      if (   idoubling == (10*IFLAG_670_220 + 0) &
        .or. idoubling == (10*IFLAG_670_220 + 1) ) then
        idoubling = IFLAG_670_220 ! separate around 410 between 220-670
      endif
      if (idoubling == iflag) then
        vph(:,:,:,ispec) = vpv(:,:,:,ispec)
        vsh(:,:,:,ispec) = vsv(:,:,:,ispec)
        eta(:,:,:,ispec) = 1.0
      endif
    enddo

    ! write new models (only vph, vsh, eta)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', vph)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', eta)

  enddo

  !====== Finalize
  if (myrank == 0) close(IOUT)

  call synchronize_all()
  call finalize_mpi()

end program
