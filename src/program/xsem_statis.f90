subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_statis "
  print '(a)', "    -  make statistical CDF/PDF of SEM gll models"
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

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_statis
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

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
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  real(dp), allocatable :: gll_volume(:,:,:,:)
  ! model
  real(dp), allocatable :: model(:,:,:,:)
  ! PDF
  real(dp) :: zmin, zmax, min_i, max_i, bin_size
  real(dp), allocatable :: z(:), pdf(:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_math: check your input arguments."
    stop
  endif

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
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! initialize arrays
  allocate(model(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll_volume(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(z(0:nbin), pdf(1:nbin))

  ! get value range
  zmax = -1.0 * huge(1.0_dp)
  zmin = huge(1.0_dp)
  do iproc = 0, (nproc-1)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    max_i = maxval(abs(model))
    min_i = minval(abs(model))
    if (max_i > zmax) zmax = max_i
    if (min_i < zmin) zmin = min_i
  enddo ! iproc

  ! compute PDF
  bin_size = (zmax - zmin) / nbin
  z = (/(i*bin_size, i=0,nbin)/)
  pdf = 0.0_dp
  do iproc = 0, (nproc-1)
    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    ! get gll_volume
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_mesh_gll_volume(mesh_data, gll_volume)
    ! collect model value into bins
    do i = 1, nbin
      pdf(i) = pdf(i) + sum(gll_volume, &
        mask=(abs(model)>=z(i-1) .and. abs(model)<z(i)) )
    enddo
  enddo
  ! normalize
  pdf = pdf / sum(pdf)

  ! print out results
  do i = 1, nbin
    !write(*, '(E12.4, E12.4, E12.4)') z(i-1), z(i), pdf(i)
    print '(E12.4, E12.4, E12.4)', z(i-1), z(i), pdf(i)
  enddo

end program
