subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_lamda_mu_to_tiso"
  print '(a)', "    - add isotropic model perturbation direction(dlamda, dmu, drho)"
  print '(a)', "      into VTI model files (vpv,vph,vsv,vsh,eta,rho)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_dmodel_lamda_mu_to_tiso \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <dmodel_dir> <out_dir>"
  print '(a)', "    <max_dlnv_allowed> <force_max_dlnv_allow> <fix_rho>"
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
  print '(a)', "  (string) out_dir:  output directory for new model"
  print '(a)', "  (float) max_dlnv_allow:  maximum velocity perturbation ratio allowed"
  print '(a)', "  (int) force_max_dlnv_allow:  flag whether max velocity perturbation should be enforced"
  print '(a)', "  (int) fix_rho:  flag whether rho is fixed"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  max_dlnv_allowed > 0: model + dmodel * step_length"
  print '(a)', "  max_dlnv_allowed < 0: model - dmodel * step_length"
  print '(a)', "  the step_length is determined when, "
  print '(a)', "  flag_force_max_dlnv = 1: force max velocity perturbation = max_dlnv_allowed"
  print '(a)', "  flag_force_max_dlnv = 0: only reduce max velocity perturbation to max_dlnv_allowed"
  print '(a)', "    if already smaller than max_dlnv_allowed, then do nothing"
  print '(a)', "  write(append) to log_file: step_length: "
end subroutine


program xsem_add_dmodel_lamda_mu_to_tiso

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir, dmodel_dir
  character(len=MAX_STRING_LEN) :: out_dir
  real(dp) :: max_dlnv_allow
  logical :: force_max_dlnv_allow, fix_rho

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vpv, vph, vsv, vsh, eta, rho
  real(dp), dimension(:,:,:,:), allocatable :: A, C, F, L, N
  ! dmodel 
  real(dp), dimension(:,:,:,:), allocatable :: dmu, dlamda, drho
  ! model perturbations
  real(dp) :: max_dln_vpv, max_dln_vph, max_dln_vsv, max_dln_vsh, max_dlnv_all
  real(dp), allocatable :: max_dlnv_procs(:)
  real(dp) :: step_length

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_dmodel_lamda_mu_to_tiso: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') model_dir 
  read(args(4),'(a)') dmodel_dir 
  read(args(5),'(a)') out_dir
  read(args(6),*) max_dlnv_allow
  select case (args(7))
    case ('0')
      force_max_dlnv_allow = .false.
    case ('1')
      force_max_dlnv_allow = .true.
    case default
      print *, '[ERROR]: force_max_dlnv_allow must be 0 or 1'
      stop
  end select
  select case (args(8))
    case ('0')
      fix_rho = .false.
    case ('1')
      fix_rho = .true.
    case default
      print *, '[ERROR]: fix_rho must be 0 or 1'
      stop
  end select

  !===== Get step length first
  allocate(max_dlnv_procs(0:nproc-1))
  max_dlnv_procs = 0.0

  ! get mesh info: nspec
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! intialize  arrays 
  allocate(vpv(NGLLX,NGLLY,NGLLZ,NSPEC), &
           vph(NGLLX,NGLLY,NGLLZ,NSPEC), &
           vsv(NGLLX,NGLLY,NGLLZ,NSPEC), &
           vsh(NGLLX,NGLLY,NGLLZ,NSPEC), &
           eta(NGLLX,NGLLY,NGLLZ,NSPEC), &
           rho(NGLLX,NGLLY,NGLLZ,NSPEC), &
          A(NGLLX,NGLLY,NGLLZ,NSPEC), &
          C(NGLLX,NGLLY,NGLLZ,NSPEC), &
          F(NGLLX,NGLLY,NGLLZ,NSPEC), &
          L(NGLLX,NGLLY,NGLLZ,NSPEC), &
          N(NGLLX,NGLLY,NGLLZ,NSPEC), &
           dmu(NGLLX,NGLLY,NGLLZ,NSPEC), &
        dlamda(NGLLX,NGLLY,NGLLZ,NSPEC), &
          drho(NGLLX,NGLLY,NGLLZ,NSPEC))
 
  ! get step length
  print *, '#====== get step length'
  do iproc = 0, (nproc-1)

    print *, '#-- iproc=', iproc

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
    max_dln_vpv = 0.5 * maxval(abs((dlamda+2.0*dmu)/C - drho/rho))
    max_dln_vph = 0.5 * maxval(abs((dlamda+2.0*dmu)/A - drho/rho))
    max_dln_vsv = 0.5 * maxval(abs(dmu/L - drho/rho))
    max_dln_vsh = 0.5 * maxval(abs(dmu/N - drho/rho))
    max_dlnv_procs(iproc) = max(max_dln_vpv,max_dln_vph,max_dln_vsv,max_dln_vsh)

    print *, "#max_dln[vpv,vph,vsv,vsh,rho]"
    print *, max_dln_vpv
    print *, max_dln_vph
    print *, max_dln_vsv
    print *, max_dln_vsh
    print *, maxval(abs(drho/rho))

  enddo ! do iproc

  ! determine the model update step length
  max_dlnv_all = maxval(max_dlnv_procs)
  if ( (max_dlnv_all>abs(max_dlnv_allow)) .or. force_max_dlnv_allow ) then
    step_length = max_dlnv_allow / max_dlnv_all
  else ! only set positive/negative sign
    step_length = SIGN(1.d0, max_dlnv_allow)
  endif

  print *, "max_dlnv_all_proc: ", max_dlnv_all
  print *, "step_length: ", step_length

  !====== make new model
  print *, "#====== make new model "
  do iproc = 0, (nproc-1)

    print *, '#-- iproc=', iproc

    ! intialize mode/dmodel arrays 
    if (.not. allocated(vpv)) then
      !deallocate(vpv,vph,vsv,vsh,eta,rho,A,C,F,L,N,dlamda,dmu,drho)
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
    endif

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

    ! calculate new models
    dlamda = dlamda * step_length 
    dmu = dmu * step_length
    drho = drho * step_length

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

end program
