subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_divide_hess_water_level "
  print '(a)', "    -  divide kernel by estimated hessian diagonals with water level"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_divide_hess_water_level \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <kernel_tags> <hess_dir> <hess_name> "
  print '(a)', "    <eps> <out_dir> <suffix>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int)    nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_<kernel_names>.bin"
  print '(a)', "  (string) kernel_tags:  parameter names, e.g. vsv2_kernel,vsh2_kernel"
  print '(a)', "  (string) hess_dir:  directory holds proc*_reg1_<hess_name>.bin"
  print '(a)', "  (string) hess_name:  name for diagonal hessian files"
  print '(a)', "  (float)  eps: water level, must between 0 and 1"
  print '(a)', "  (string) out_dir:  output directory "
  print '(a)', "  (string) suffix:   proc*_reg1_<kernel_names>_<suffix>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. out = kernel / (diag_hess_normalized + eps)"
  print '(a)', "  3. diag_hess is normalized to maximum amplitude of 1, eps can be chosen by trial and error (e.g. ~ 0.001)"
  print '(a)', "  4. it's better to smooth the hessian before use in order to have smoothed hess inverse"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_kernel_divide_hess_water_level
  
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
  character(len=MAX_STRING_LEN) :: kernel_dir, kernel_tags
  character(len=MAX_STRING_LEN) :: hess_dir, hess_name
  real(dp) :: eps 
  character(len=MAX_STRING_LEN) :: out_dir, out_suffix

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! kernel 
  integer :: iker, nkernel
  character(len=MAX_STRING_LEN), allocatable :: kernel_names(:)
  real(dp), allocatable :: kernel_gll(:,:,:,:)
  ! hess
  real(dp) :: max_hess_local, max_hess_all
  real(dp), allocatable :: hess_gll(:,:,:,:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_kernel_divide_hess_water_level: check your input arguments."
      call abort_mpi()
    endif 
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') kernel_dir
  read(args(4), '(a)') kernel_tags
  read(args(5), '(a)') hess_dir
  read(args(6), '(a)') hess_name
  read(args(7), *) eps
  read(args(8), '(a)') out_dir
  read(args(9), '(a)') out_suffix

  ! validate inputs
  if (myrank == 0) then
    if (eps <= 0.d0 .or. eps >= 1.0d0) then
      print *, "[ERROR] eps must beween 0 and 1."
      call abort_mpi()
    endif
  endif 
  print *, "eps = ", eps

  !===== loop each mesh/model slice

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! get number of kernel names
  call sem_utils_delimit_string(kernel_tags, ',', kernel_names, nkernel)

  ! initialize arrays
  allocate(kernel_gll(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(hess_gll(NGLLX,NGLLY,NGLLZ,nspec))

  ! get max value of hess
  max_hess_local = 0.0_dp
  do iproc = myrank, (nproc-1), nrank
    call sem_io_read_gll_file_1(hess_dir, iproc, iregion, hess_name, hess_gll)
    max_hess_local = max(max_hess_local, maxval(abs(hess_gll)))
    print *, "max_hess_local = ", max_hess_local
  enddo
  call synchronize_all()
  call max_all_dp(max_hess_local, max_hess_all)
  call bcast_all_singledp(max_hess_all)
  call synchronize_all()

  print *, "max_hess_all = ", max_hess_all

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc = ', iproc

    ! read hess
    call sem_io_read_gll_file_1(hess_dir, iproc, iregion, hess_name, hess_gll)
    ! normalize hess
    hess_gll = hess_gll / max_hess_all

    do iker = 1, nkernel
      ! read kernel
      call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, kernel_names(iker), kernel_gll)
      ! water level division 
      kernel_gll = kernel_gll / (hess_gll + eps)
      ! write out
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, trim(kernel_names(iker))//"_"//out_suffix, kernel_gll)
    enddo

  enddo

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
