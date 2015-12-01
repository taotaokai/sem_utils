subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_lamda_mu_rho_to_tiso"
  print '(a)', "    - add isotropic model perturbation direction(dlamda, dmu, drho)"
  print '(a)', "      into VTI model files (vpv,vph,vsv,vsh,eta,rho)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_lamda_mu_to_tiso \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <dmodel_dir> "
  print '(a)', "    <max_dlnv_allowed> <force_max_dlnv_allow> <fix_rho>"
  print '(a)', "    <out_dir> <log_file>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  new_model = model + step_length * dmodel"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
  print '(a)', "  (string) dmodel_dir:  directory holds proc*_reg1_[lamda_dmodel,mu_dmodel,rho_dmodel].bin"
  print '(a)', "  (float) max_dlnv_allow:  maximum velocity perturbation ratio allowed"
  print '(a)', "  (int) force_max_dlnv_allow:  flag whether max velocity perturbation should be enforced"
  print '(a)', "  (int) fix_rho:  flag whether rho is fixed"
  print '(a)', "  (string) out_dir:  output directory for new model"
  print '(a)', "  (string) log_file:  file to log runing info"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. max_dlnv_allowed > 0: model + dmodel * step_length"
  print '(a)', "     max_dlnv_allowed < 0: model - dmodel * step_length"
  print '(a)', "  2. the step_length is determined when, "
  print '(a)', "    flag_force_max_dlnv = 1: scale to have max velocity perturbation = max_dlnv_allowed"
  print '(a)', "    flag_force_max_dlnv = 0: only effective when max velocity perturbation is larger than max_dlnv_allowed"
  print '(a)', "  3. can be run in parallel"

end subroutine


program xsem_add_dmodel_lamda_mu_to_tiso

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
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: dmodel_dir
  real(dp) :: max_dlnv_allow
  logical :: force_max_dlnv_allow
  logical :: fix_rho
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: log_file

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vpv, vph, vsv, vsh, eta, rho
  real(dp), dimension(:,:,:,:), allocatable :: A, C, F, L, N
  ! dmodel 
  real(dp), dimension(:,:,:,:), allocatable :: dmu, dlamda, drho
  ! model perturbations
  real(dp) :: max_dln_vpv, max_dln_vph, max_dln_vsv, max_dln_vsh
  real(dp) :: max_dln_vpv_all, max_dln_vph_all, max_dln_vsv_all, max_dln_vsh_all
  real(dp) :: max_dlnv_all
  real(dp) :: step_length

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_add_dmodel_lamda_mu_rho_to_tiso: check your inputs."
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
  read(args(4),'(a)') dmodel_dir 
  read(args(5),*) max_dlnv_allow
  select case (args(6))
    case ('0')
      force_max_dlnv_allow = .false.
    case ('1')
      force_max_dlnv_allow = .true.
    case default
      if (myrank==0) then
        print *, '[ERROR]: force_max_dlnv_allow must be 0 or 1'
        call abort_mpi()
      endif
  end select
  select case (args(7))
    case ('0')
      fix_rho = .false.
    case ('1')
      fix_rho = .true.
    case default
      if (myrank==0) then
        print *, '[ERROR]: fix_rho must be 0 or 1'
        call abort_mpi()
      endif
  end select
  read(args(8),'(a)') out_dir
  read(args(9),'(a)') log_file

  call synchronize_all()

  !===== Get step length first

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(vpv(NGLLX,NGLLY,NGLLZ,nspec), &
           vph(NGLLX,NGLLY,NGLLZ,nspec), &
           vsv(NGLLX,NGLLY,NGLLZ,nspec), &
           vsh(NGLLX,NGLLY,NGLLZ,nspec), &
           eta(NGLLX,NGLLY,NGLLZ,nspec), &
           rho(NGLLX,NGLLY,NGLLZ,nspec), &
             A(NGLLX,NGLLY,NGLLZ,nspec), &
             C(NGLLX,NGLLY,NGLLZ,nspec), &
             F(NGLLX,NGLLY,NGLLZ,nspec), &
             L(NGLLX,NGLLY,NGLLZ,nspec), &
             N(NGLLX,NGLLY,NGLLZ,nspec), &
           dmu(NGLLX,NGLLY,NGLLZ,nspec), &
        dlamda(NGLLX,NGLLY,NGLLZ,nspec), &
          drho(NGLLX,NGLLY,NGLLZ,nspec))
 
  ! get maximum velocity perturbation from dmodel
  max_dln_vpv = 0.0_dp
  max_dln_vph = 0.0_dp
  max_dln_vsv = 0.0_dp
  max_dln_vsh = 0.0_dp

  ! open log file for write
  if (myrank == 0) then

    open(IOUT, file=trim(log_file), status='unknown', &
      form='formatted', action='write', iostat=ier)

    if (ier /= 0) then
      write(*,*) '[ERROR] xsem_add_dmodel_lamda_mu_to_tiso: failed to open file ', trim(log_file)
      call abort_mpi()
    endif

  endif

  call synchronize_all()

  if (myrank ==0) print '(a)', "#====== get step_length"

  do iproc = myrank, (nproc-1), nrank

    !if (myrank ==0) write(IOUT,'(a,2X,I4)') "# iproc=", iproc

    ! read old models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    ! get A,C,F,L,N from old model
    A = vph**2 * rho
    C = vpv**2 * rho
    L = vsv**2 * rho
    N = vsh**2 * rho
    F = eta * (A - 2.0*L)
    ! read dmodel: update direction 
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'lamda_dmodel', dlamda)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'mu_dmodel', dmu)
    if (fix_rho) then
      drho = 0.0_dp 
    else
      call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'rho_dmodel', drho)
    endif
    ! calculate maximum absolute relative perturbation of velocities
    max_dln_vpv = max(max_dln_vpv, 0.5*maxval(abs((dlamda+2.0*dmu)/C - drho/rho)))
    max_dln_vph = max(max_dln_vph, 0.5*maxval(abs((dlamda+2.0*dmu)/A - drho/rho)))
    max_dln_vsv = max(max_dln_vsv, 0.5*maxval(abs(dmu/L - drho/rho)))
    max_dln_vsh = max(max_dln_vsh, 0.5*maxval(abs(dmu/N - drho/rho)))
  enddo ! do iproc

  call synchronize_all()

  call max_all_dp(max_dln_vpv, max_dln_vpv_all)
  call max_all_dp(max_dln_vph, max_dln_vph_all)
  call max_all_dp(max_dln_vsv, max_dln_vsv_all)
  call max_all_dp(max_dln_vsh, max_dln_vsh_all)

  ! determine step_length
  if (myrank == 0) then

    max_dlnv_all = max(max_dln_vpv_all,max_dln_vph_all, &
                       max_dln_vsv_all,max_dln_vsh_all)

    if (max_dlnv_all>abs(max_dlnv_allow) .or. force_max_dlnv_allow) then
      step_length = max_dlnv_allow / max_dlnv_all
    else 
      ! only set positive/negative sign
      step_length = SIGN(1.d0, max_dlnv_allow)
    endif

    write(IOUT,'(a)') "#[LOG] xsem_add_dmodel_lamda_mu_to_tiso"
    write(IOUT,'(a,2X,E12.4)') "max_dln_vpv_all=", max_dln_vpv_all
    write(IOUT,'(a,2X,E12.4)') "max_dln_vph_all=", max_dln_vph_all
    write(IOUT,'(a,2X,E12.4)') "max_dln_vsv_all=", max_dln_vsv_all
    write(IOUT,'(a,2X,E12.4)') "max_dln_vsh_all=", max_dln_vsh_all
    write(IOUT,'(a,2X,E12.4)') "max_dlnv_all=", max_dlnv_all
    write(IOUT,'(a,2X,E12.4)') "step_length=", step_length

    close(IOUT)

  endif

  call synchronize_all()

  call bcast_all_singledp(step_length)

  !====== create new model
  if (myrank == 0) print '(a)', "#====== create new model"

  do iproc = myrank, (nproc-1), nrank

    print '(a,2X,I4)', "# iproc=", iproc

    ! read old models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', eta)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    ! get A,C,F,L,N from old model
    A = vph**2 * rho; 
    C = vpv**2 * rho;
    L = vsv**2 * rho; 
    N = vsh**2 * rho;
    F = eta * (A - 2.0*L)
    ! read dmodel: update direction 
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'lamda_dmodel', dlamda)
    call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'mu_dmodel', dmu)
    if (fix_rho) then
      drho = 0.0_dp
    else
      call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'rho_dmodel', drho)
    endif
    ! scale dmodel by step_length
    dlamda = dlamda * step_length 
    dmu = dmu * step_length
    drho = drho * step_length
    ! add dmodel
    rho = rho + drho
    vpv = sqrt((C+dlamda+2.0*dmu)/rho)
    vph = sqrt((A+dlamda+2.0*dmu)/rho)
    vsv = sqrt((L+dmu)/rho)
    vsh = sqrt((N+dmu)/rho)
    eta = (F+dlamda)/(A-2.0*L+dlamda)
    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', vpv)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', vph)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', vsv)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', vsh)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', eta)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho', rho)

  enddo

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
