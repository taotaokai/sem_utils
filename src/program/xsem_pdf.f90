subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_pdf "
  print '(a)', "    -  compute PDF of volumetric integral of model amplitude"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_statis \"
  print '(a)', "    <nproc> <mesh_dir>"
  print '(a)', "    <model_dir> <model_name> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_name>.bin"
  print '(a)', "  (string) model_name:  model name, e.g. mu_kernel "
  print '(a)', "  (int) nbin:  number of bins (from min to max value)"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. output is a list of <bin_left> <bin_right> <pdf>"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_pdf
  
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
  character(len=MAX_STRING_LEN) :: model_dir, model_name
  integer :: nbin

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  real(dp), allocatable :: gll_volume(:,:,:,:)
  ! model
  real(dp), allocatable :: model(:,:,:,:)
  ! PDF
  real(dp) :: zmin, zmax, min_i, max_i, bin_size
  real(dp) :: zmin_all, zmax_all
  real(dp), allocatable :: z(:), pdf(:), pdf_all(:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_pdf: check your input arguments."
      call abort_mpi()
    endif 
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_name
  read(args(5), *) nbin

  !===== loop each mesh/model slice

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize arrays
  allocate(model(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll_volume(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(z(0:nbin), pdf(nbin))
  if (myrank == 0) allocate(pdf_all(nbin))

  ! get model value range
  zmax = -1.0 * huge(1.0_dp)
  zmin = huge(1.0_dp)
  do iproc = myrank, (nproc-1), nrank
    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    max_i = maxval(model)
    min_i = minval(model)
    if (max_i > zmax) zmax = max_i
    if (min_i < zmin) zmin = min_i
  enddo ! iproc

  call synchronize_all()
   
  call max_all_dp(zmax, zmax_all)
  call min_all_dp(zmin, zmin_all)
  call bcast_all_singledp(zmax_all)
  call bcast_all_singledp(zmin_all)

  call synchronize_all()

  if (myrank == 0) print *, '# min/max=', zmin_all, zmax_all

  ! compute PDF
  bin_size = (zmax_all - zmin_all) / nbin
  z = (/(zmin_all + i*bin_size, i=0,nbin)/)
  pdf = 0.0_dp

  do iproc = myrank, (nproc-1), nrank
    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    ! get gll_volume
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_mesh_gll_volume(mesh_data, gll_volume)
    ! get volumetric integral of model amplitudes within each bin 
    do i = 1, nbin
      pdf(i) = pdf(i) + sum(gll_volume*abs(model), &
        mask=(model>=z(i-1) .and. model<z(i)) )
    enddo
  enddo
  
  call synchronize_all()

  call sum_all_dp_n(pdf, pdf_all, nbin)

  call synchronize_all()

  ! normalize PDF
  if (myrank == 0) then
    print *, "# Integral(model,dV)=", sum(pdf_all)
    pdf_all = pdf_all / sum(pdf_all)
  endif

  ! print out results
  if (myrank == 0) then
    do i = 1, nbin
      print '(E12.4, E12.4, E12.4)', z(i-1), z(i), pdf_all(i)
    enddo
  endif

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program