subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_pcg_dmodel"
  print '(a)', "    - calculate precondtioned conjugate gradient update direction"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_pcg_dmodel"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <current_kernel_dir> <current_kernel_name> "
  print '(a)', "    <previous_kernel_dir> <previous_kernel_name> "
  print '(a)', "    <previous_dmodel_dir> <previous_dmodel_name> "
  print '(a)', "    <out_dir> <out_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  To maximize an objective function the PCG update is d_k+1 = p_k+1 - beta_k*d_k"
  print '(a)', "  , where beta_k = < y_k - 2*<y_k,y_k>/<d_k,y_k>d_k , p_k+1 > / <d_k,y_k>"
  print '(a)', "  (Hager and Zhang, 2003)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) current_kernel_dir: kernel directory for the current model"
  print '(a)', "  (string) current_kernel_name:  proc*_reg1_<kernel_name>.bin"
  print '(a)', "  (string) previous_kernel_dir: kernel directory for the previous model"
  print '(a)', "  (string) previous_kernel_name:  proc*_reg1_<kernel_name>.bin"
  print '(a)', "  (string) previous_dmodel_dir: dmodel directory for the previous model"
  print '(a)', "  (string) previous_dmodel_name:  proc*_reg1_<dmodel_name>.bin"
  print '(a)', "  (string) out_dir:  output directory for dmodel files"
  print '(a)', "  (string) out_name:  proc*_reg1_<out_name>.bin"
  print '(a)', ""
  print '(a)', "note"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

END SUBROUTINE


program xsem_kernel_vti_3pars_pcg_dmodel

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 10
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: current_kernel_dir
  character(len=MAX_STRING_LEN) :: current_kernel_name
  character(len=MAX_STRING_LEN) :: previous_kernel_dir
  character(len=MAX_STRING_LEN) :: previous_kernel_name
  character(len=MAX_STRING_LEN) :: previous_dmodel_dir
  character(len=MAX_STRING_LEN) :: previous_dmodel_name
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! kernel
  real(dp), dimension(:,:,:,:), allocatable :: kernel_current
  real(dp), dimension(:,:,:,:), allocatable :: kernel_previous
  ! previous dmodel
  real(dp), dimension(:,:,:,:), allocatable :: dmodel
  ! cg update parameter
  real(dp), dimension(:,:,:,:), allocatable :: dkernel
  real(dp), dimension(:,:,:,:), allocatable :: volume_gll
  real(dp) :: yk_dk,  yk_dk_all
  real(dp) :: yk_pk1, yk_pk1_all
  real(dp) :: yk_yk,  yk_yk_all
  real(dp) :: dk_pk1, dk_pk1_all
  real(dp) :: beta

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
  read(args(3),'(a)') current_kernel_dir 
  read(args(4),'(a)') current_kernel_name
  read(args(5),'(a)') previous_kernel_dir 
  read(args(6),'(a)') previous_kernel_name
  read(args(7),'(a)') previous_dmodel_dir 
  read(args(8),'(a)') previous_dmodel_name
  read(args(9),'(a)') out_dir
  read(args(10),'(a)') out_name

  call synchronize_all()

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays
  allocate(kernel_current(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kernel_previous(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dmodel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dkernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(volume_gll(NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate inner products
  yk_dk = 0.0_dp
  yk_pk1 = 0.0_dp
  yk_yk = 0.0_dp
  dk_pk1 = 0.0_dp

  do iproc = myrank, (nproc-1), nrank
  
    ! read in mesh 
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    ! calculate gll volumes
    call sem_mesh_gll_volume(mesh_data, volume_gll)

    ! read kernels (p_k+1)
    call sem_io_read_gll_file_1(current_kernel_dir, iproc, iregion,  current_kernel_name, kernel_current)

    ! read kernels (p_k)
    call sem_io_read_gll_file_1(previous_kernel_dir, iproc, iregion, previous_kernel_name, kernel_previous)

    ! read previous dmodel (d_k)
    call sem_io_read_gll_file_1(previous_dmodel_dir, iproc, iregion, previous_dmodel_name, dmodel)

    ! y_k = p_k+1 - p_k
    dkernel = kernel_current - kernel_previous

    ! <y_k, d_k>
    yk_dk = yk_dk + sum(dkernel*dmodel*volume_gll)

    ! <y_k, p_k+1>
    yk_pk1 = yk_pk1 + sum(dkernel*kernel_current*volume_gll)

    ! <y_k, y_k>
    yk_yk = yk_yk + sum(dkernel**2 * volume_gll)

    ! <d_k, p_k+1>
    dk_pk1 = dk_pk1 + sum(dmodel*kernel_current*volume_gll)

  enddo

  call synchronize_all()

  call sum_all_dp(yk_dk, yk_dk_all)
  call sum_all_dp(yk_pk1, yk_pk1_all)
  call sum_all_dp(yk_yk, yk_yk_all)
  call sum_all_dp(dk_pk1, dk_pk1_all)

  call bcast_all_singledp(yk_dk_all)
  call bcast_all_singledp(yk_pk1_all)
  call bcast_all_singledp(yk_yk_all)
  call bcast_all_singledp(dk_pk1_all)

  print *, "yk_dk = ", yk_dk_all
  print *, "yk_pk1 = ", yk_pk1_all
  print *, "yk_yk = ", yk_yk_all
  print *, "dk_pk1 = ", dk_pk1_all

  !====== make new dmodel

  call synchronize_all()

  ! cg update parameter
  !beta = (yk_pk1_all - 2.0*dk_pk1_all*yk_yk_all/yk_dk_all) / yk_dk_all ! Hager and Zhang (2003)
  beta = yk_pk1_all/yk_dk_all ! Hestenes and Stiefel (1952)

  print *, "update parameter = ", beta 

  do iproc = myrank, (nproc-1), nrank

    print *, "iproc = ", iproc

    ! read kernels (p_k+1)
    call sem_io_read_gll_file_1(current_kernel_dir, iproc, iregion,  current_kernel_name, kernel_current)

    ! read previous dmodel (d_k)
    call sem_io_read_gll_file_1(previous_dmodel_dir, iproc, iregion, previous_dmodel_name, dmodel)

    ! get new dmodel
    dmodel = kernel_current - beta*dmodel

    ! write new dmodel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, out_name, dmodel)

  enddo
 
  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
