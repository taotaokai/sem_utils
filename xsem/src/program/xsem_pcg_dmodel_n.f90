subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_pcg_dmodel_n"
  print '(a)', "    - calculate (precondtioned) conjugate gradient update direction"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_pcg_dmodel"
  print '(a)', "    <nproc> <mesh_dir> "
  print '(a)', "    <current_precond_kernel_dir> <current_precond_kernel_names> "
  print '(a)', "    <current_kernel_dir> <current_kernel_names> "
  print '(a)', "    <previous_kernel_dir> <previous_kernel_names> "
  print '(a)', "    <previous_dmodel_dir> <previous_dmodel_names> "
  print '(a)', "    <cg_type> <out_dir> <out_names>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  To maximize an objective function the PCG update is d_k+1 = p_k+1 - beta_k*d_k"
  print '(a)', "  , where beta_k = < y_k - 2*<y_k,y_k>/<d_k,y_k>d_k , p_k+1 > / <d_k,y_k>, y_k = g_k+1 - g_k"
  print '(a)', "  (Hager and Zhang, 2003)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) current_precond_kernel_dir: directory of the current preconditioned kernel files proc*_reg1_<kernel_names>.bin"
  print '(a)', "  (string) current_precond_kernel_names:  comma separated string, e.g. dlnvs_kernel,kappa_kernel,... "
  print '(a)', "  (string) previous_kernel_dir: directory for the previous kernel files"
  print '(a)', "  (string) current_kernel_dir: directory of the current kernel files proc*_reg1_<kernel_names>.bin"
  print '(a)', "  (string) current_kernel_names:  comma separated string, e.g. dlnvs_kernel,kappa_kernel,... "
  print '(a)', "  (string) previous_kernel_dir: directory for the previous kernel files"
  print '(a)', "  (string) previous_kernel_names:  comma separated string"
  print '(a)', "  (string) previous_dmodel_dir: directory for the previous dmodel (model update direction)"
  print '(a)', "  (string) previous_dmodel_names:  comma separated string"
  print '(a)', "  (string) cg_type:  choise on update parameter, e.g. HS, N"
  print '(a)', "  (string) out_dir:  output directory for new dmodel files"
  print '(a)', "  (string) out_names:  comma separated string"
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
  integer, parameter :: nargs = 13
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: current_precond_kernel_dir
  character(len=MAX_STRING_LEN) :: current_precond_kernel_names
  character(len=MAX_STRING_LEN) :: current_kernel_dir
  character(len=MAX_STRING_LEN) :: current_kernel_names
  character(len=MAX_STRING_LEN) :: previous_kernel_dir
  character(len=MAX_STRING_LEN) :: previous_kernel_names
  character(len=MAX_STRING_LEN) :: previous_dmodel_dir
  character(len=MAX_STRING_LEN) :: previous_dmodel_names
  character(len=MAX_STRING_LEN) :: cg_type 
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_names

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model names
  integer :: nmodel, nmodel_1, nmodel_2, nmodel_3, nmodel_4
  character(len=MAX_STRING_LEN), allocatable :: current_precond_kernel_name_list(:)
  character(len=MAX_STRING_LEN), allocatable :: current_kernel_name_list(:)
  character(len=MAX_STRING_LEN), allocatable :: previous_kernel_name_list(:)
  character(len=MAX_STRING_LEN), allocatable :: previous_dmodel_name_list(:)
  character(len=MAX_STRING_LEN), allocatable :: out_name_list(:)
  ! kernel
  real(dp), dimension(:,:,:,:,:), allocatable :: pk1 !current preconditioned kernel
  real(dp), dimension(:,:,:,:,:), allocatable :: gk1 !current gradient/kernel
  real(dp), dimension(:,:,:,:,:), allocatable :: gk  !previous gradient
  ! previous dmodel
  real(dp), dimension(:,:,:,:,:), allocatable :: dk  !previous model update direction
  ! cg update parameter
  real(dp), dimension(:,:,:,:,:), allocatable :: yk  !gradient difference
  real(dp), dimension(:,:,:,:), allocatable :: volume_gll
  real(dp) :: yk_dk,  yk_dk_all
  real(dp) :: yk_pk1, yk_pk1_all
  real(dp) :: yk_yk,  yk_yk_all
  real(dp) :: dk_pk1, dk_pk1_all
  real(dp) :: dk_gk1, dk_gk1_all
  real(dp) :: dk_dk, dk_dk_all
  real(dp) :: gk1_gk1, gk1_gk1_all
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
  read(args(3),'(a)') current_precond_kernel_dir 
  read(args(4),'(a)') current_precond_kernel_names
  read(args(5),'(a)') current_kernel_dir 
  read(args(6),'(a)') current_kernel_names
  read(args(7),'(a)') previous_kernel_dir 
  read(args(8),'(a)') previous_kernel_names
  read(args(9),'(a)') previous_dmodel_dir 
  read(args(10),'(a)') previous_dmodel_names
  read(args(11),'(a)') cg_type
  read(args(12),'(a)') out_dir
  read(args(13),'(a)') out_names

  call synchronize_all()

  ! read in model names
  call sem_utils_delimit_string(current_precond_kernel_names, ',', current_precond_kernel_name_list, nmodel)
  call sem_utils_delimit_string(current_kernel_names, ',', current_kernel_name_list, nmodel_1)
  call sem_utils_delimit_string(previous_kernel_names, ',', previous_kernel_name_list, nmodel_2)
  call sem_utils_delimit_string(previous_dmodel_names, ',', previous_dmodel_name_list, nmodel_3)
  call sem_utils_delimit_string(out_names, ',', out_name_list, nmodel_4)
  if (nmodel_1 /= nmodel .or. nmodel_2 /= nmodel .or. nmodel_3 /= nmodel .or. nmodel_4 /= nmodel) then
    if (myrank == 0) then
      print *, "[ERROR] number of parameters should be the same"
      call abort_mpi()
    endif
  endif

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays
  allocate(volume_gll(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(pk1(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gk1(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gk(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dk(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
  allocate(yk(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

  !====== calculate inner products
  yk_dk = 0.0_dp
  yk_pk1 = 0.0_dp
  yk_yk = 0.0_dp
  dk_pk1 = 0.0_dp

  dk_gk1 = 0.0_dp
  dk_dk = 0.0_dp
  gk1_gk1 = 0.0_dp

  do iproc = myrank, (nproc-1), nrank
  
    ! read in mesh 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    ! calculate gll volumes
    call sem_mesh_gll_volume(mesh_data, volume_gll)

    ! read preconditioned kernel (p_k+1)
    call sem_io_read_gll_file_n(current_precond_kernel_dir, iproc, iregion, current_precond_kernel_name_list, nmodel, pk1)

    ! read kernels (g_k+1)
    call sem_io_read_gll_file_n(current_kernel_dir, iproc, iregion, current_kernel_name_list, nmodel, gk1)

    ! read kernels (g_k)
    call sem_io_read_gll_file_n(previous_kernel_dir, iproc, iregion, previous_kernel_name_list, nmodel, gk)

    ! read previous dmodel (d_k)
    call sem_io_read_gll_file_n(previous_dmodel_dir, iproc, iregion, previous_dmodel_name_list, nmodel, dk)

    ! y_k = g_k+1 - g_k
    yk = gk1 - gk

    ! <y_k, d_k>
    yk_dk = yk_dk + sum(yk*dk*spread(volume_gll,1,nmodel))

    ! <y_k, p_k+1>
    yk_pk1 = yk_pk1 + sum(yk*pk1*spread(volume_gll,1,nmodel))

    ! <y_k, y_k>
    yk_yk = yk_yk + sum(yk**2 * spread(volume_gll,1,nmodel))

    ! <d_k, p_k+1>
    dk_pk1 = dk_pk1 + sum(dk*pk1*spread(volume_gll,1,nmodel))

    ! <d_k, g_k+1>
    dk_gk1 = dk_gk1 + sum(dk*gk1 * spread(volume_gll,1,nmodel))

    ! <d_k, d_k>
    dk_dk = dk_dk + sum(dk**2 * spread(volume_gll,1,nmodel))

    ! <g_k+1, g_k+1>
    gk1_gk1 = gk1_gk1 + sum(gk1**2 * spread(volume_gll,1,nmodel))

  enddo

  call synchronize_all()

  call sum_all_dp(yk_dk, yk_dk_all)
  call sum_all_dp(yk_pk1, yk_pk1_all)
  call sum_all_dp(yk_yk, yk_yk_all)
  call sum_all_dp(dk_pk1, dk_pk1_all)

  call sum_all_dp(dk_gk1, dk_gk1_all)
  call sum_all_dp(dk_dk, dk_dk_all)
  call sum_all_dp(gk1_gk1, gk1_gk1_all)

  call bcast_all_singledp(yk_dk_all)
  call bcast_all_singledp(yk_pk1_all)
  call bcast_all_singledp(yk_yk_all)
  call bcast_all_singledp(dk_pk1_all)

  if (myrank == 0) then
    print *, "<y_k, d_k> = ", yk_dk_all
    print *, "<y_k, p_k+1> = ", yk_pk1_all
    print *, "|y_k|^2 = ", yk_yk_all
    print *, "<d_k, p_k+1> = ", dk_pk1_all

    print *, "<d_k, g_k+1> = ", dk_gk1_all
    print *, "|d_k|^2 = ", dk_dk_all
    print *, "|g_k+1|^2 = ", gk1_gk1_all
    print *, "<d_k, g_k+1>/(|d_k|*|g_k+1|) = ", dk_gk1_all/sqrt(dk_dk_all)/sqrt(gk1_gk1_all)
  endif

  !====== make new dmodel

  call synchronize_all()

  ! cg update parameter
  if (cg_type == "HS") then
    beta = yk_pk1_all/yk_dk_all ! Hestenes and Stiefel (1952)
  !else if (cg_type == "N") then
  !  beta = (yk_pk1_all - 2.0*dk_pk1_all*yk_yk_all/yk_dk_all) / yk_dk_all ! Hager and Zhang (2003)
  else
    print *, "[ERROR] unknown type of CG (valid option: HS or N)"
    call abort_mpi()
  endif

  if (myrank == 0) then
    print *, "update parameter = ", beta 
  endif

  do iproc = myrank, (nproc-1), nrank

    !print *, "iproc = ", iproc

    ! read kernels (p_k+1)
    call sem_io_read_gll_file_n(current_precond_kernel_dir, iproc, iregion, current_precond_kernel_name_list, nmodel, pk1)

    ! read previous dmodel (d_k)
    call sem_io_read_gll_file_n(previous_dmodel_dir, iproc, iregion, previous_dmodel_name_list, nmodel, dk)

    ! get new dmodel (d_k+1)
    ! note: because I am trying to maximize the objective function, the update formula is d_k = g_k - beta*d_k-1
    dk = pk1 - beta*dk

    ! write new dmodel
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, out_name_list, nmodel, dk)

  enddo
 
  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
