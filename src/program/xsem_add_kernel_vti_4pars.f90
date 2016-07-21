subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_kernel_vti_4pars"
  print '(a)', "    - add vti (precondioned-)kernel (vphi2, vsv2, vsh2, vf2)"
  print '(a)', "      into vti model file (vpv,vph,vsv,vsh,eta,rho)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_vti_kernel \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <kernel_suffix> "
  print '(a)', "    <max_dlnv_allowed> <flag_force_max_dlnv> <ratio_rho>"
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
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_[vphi2,vsv2,vsh2,vf2]_kernel.bin"
  print '(a)', "  (string) kernel_suffix:  proc*_reg1_[vphi2,vsv2,vsh2,vf2]_kernel<kernel_suffix>.bin"
  print '(a)', "  (float) max_dlnv_allow:  maximum velocity perturbation ratio allowed"
  print '(a)', "  (int) force_max_dlnv_allowed:  flag whether max velocity perturbation should be enforced"
  print '(a)', "  (float) ratio_rho:  dln(rho) = ratio_rho * dln(vs_avg)"
  print '(a)', "  (string) out_dir:  output directory for new model"
  print '(a)', "  (string) log_file:  file to log runing info"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. max_dlnv_allowed > 0: model + dmodel * step_length"
  print '(a)', "     max_dlnv_allowed < 0: model - dmodel * step_length"
  print '(a)', "  3. the step_length is determined when, "
  print '(a)', "     flag_force_max_dlnv = 1: scale to have max velocity perturbation = max_dlnv_allowed"
  print '(a)', "     flag_force_max_dlnv = 0: only effective when max velocity perturbation is larger than max_dlnv_allowed"
  print '(a)', "  4. vsv2 = vsv^2, vs_avg = (vsh+vsv)/2"

end subroutine


program xsem_add_vti_kernel

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
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: kernel_suffix
  real(dp) :: max_dlnv_allowed
  logical :: flag_force_max_dlnv
  real(dp) :: ratio_rho
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
  real(dp), dimension(:,:,:,:), allocatable :: vphi2, vpv2, vph2, vsv2, vsh2, vf2, rho
  ! kernel
  real(dp), dimension(:,:,:,:), allocatable :: vphi2_kernel, vsv2_kernel, vsh2_kernel, vf2_kernel
  ! step length
  real(dp), dimension(:,:,:,:), allocatable :: dln_rho
  real(dp) :: max_dlnv, max_dlnv_all
  real(dp) :: step_length

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_add_vti_kernel: check your inputs."
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
  read(args(4),'(a)') kernel_dir 
  read(args(5),'(a)') kernel_suffix
  read(args(6),*) max_dlnv_allowed
  select case (args(7))
    case ('0')
      flag_force_max_dlnv = .false.
    case ('1')
      flag_force_max_dlnv = .true.
    case default
      if (myrank==0) then
        print *, '[ERROR]: flag_force_max_dlnv must be 0 or 1'
        call abort_mpi()
      endif
  end select
  read(args(8),*) ratio_rho 
  read(args(9),'(a)') out_dir
  read(args(10),'(a)') log_file

  call synchronize_all()

  !====== open log file for write
  if (myrank == 0) then

    open(IOUT, file=trim(log_file), status='unknown', &
      form='formatted', action='write', iostat=ier)

    if (ier /= 0) then
      write(*,*) '[ERROR] xsem_add_dmodel_lamda_mu_to_tiso: failed to open file ', trim(log_file)
      call abort_mpi()
    endif

    write(IOUT, '(a)') '#[LOG] xsem_add_dmodel_lamda_mu_to_tiso'

  endif

  !===== Get step length first

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(vphi2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vph2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vf2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vphi2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vf2_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dln_rho(NGLLX,NGLLY,NGLLZ,nspec))
 
  ! get maximum velocity perturbation from dmodel
  max_dlnv = 0.0_dp

  if (myrank ==0) write(IOUT,'(a)') "#====== get step_length"
  call synchronize_all()

  do iproc = myrank, (nproc-1), nrank

    ! read models
    !call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph2)
    !call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh2)
    !call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', vf2)
    !call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    ! convert to squared value and also vf2
    !vph2 = vph2**2
    !vpv2 = vpv2**2
    vsv2 = vsv2**2
    vsh2 = vsh2**2
    !vf2 = vf2*(vph2 - 2.0*vsv2)

    ! read kernels
    !call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'vphi2_kernel'//kernel_suffix, vphi2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsv2_kernel'//kernel_suffix, vsv2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsh2_kernel'//kernel_suffix, vsh2_kernel)
    !call sem_io_read_gll_file_1(dmodel_dir, iproc, iregion, 'vf2_kernel'//kernel_suffix, vphi2_kernel)

    ! calculate maximum absolute relative perturbation of velocities, only use vsv and vsh
    max_dlnv = max(max_dlnv, 0.5*maxval(abs(vsv2_kernel)/vsv2), 0.5*maxval(abs(vsh2_kernel)/vsh2))
  enddo ! do iproc

  call synchronize_all()
  call max_all_dp(max_dlnv, max_dlnv_all)

  ! determine step_length
  if (myrank == 0) then
    if (max_dlnv_all>abs(max_dlnv_allowed) .or. flag_force_max_dlnv) then
      step_length = max_dlnv_allowed / max_dlnv_all
    else 
      ! only set positive/negative sign
      step_length = SIGN(1.d0, max_dlnv_allowed)
    endif
  endif

  ! log 
  if (myrank == 0) then
    write(IOUT,'(a,2X,E12.4)') "max_dlnv_all=", max_dlnv_all
    write(IOUT,'(a,2X,E12.4)') "step_length=", step_length
  endif

  call synchronize_all()
  call bcast_all_singledp(step_length)

  !====== create new model
  if (myrank == 0) write(IOUT,'(a)') "#====== create new model"
  call synchronize_all()

  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', vf2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    ! convert to squared value and also vf2
    vph2 = vph2**2
    vpv2 = vpv2**2
    vsv2 = vsv2**2
    vsh2 = vsh2**2
    vf2 = vf2*(vph2 - 2.0*vsv2)
    vphi2 = (vph2 + vpv2 - 4.0/3.0*(vsh2+vsv2))/2.0 ! to enforce Vphi,Vsv,Vsh,F parameterization of the output model 

    ! read kernels
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vphi2_kernel'//kernel_suffix, vphi2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsv2_kernel'//kernel_suffix, vsv2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsh2_kernel'//kernel_suffix, vsh2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vf2_kernel'//kernel_suffix, vf2_kernel)

    ! scale
    dln_rho = ratio_rho*step_length*0.25*(vsv2_kernel/vsv2 + vsh2_kernel/vsh2)
    rho = rho*(1.0 + dln_rho)
    print *, "dln_rho: min/max = ", minval(dln_rho), maxval(dln_rho)
    vphi2 = vphi2 + step_length*vphi2_kernel
    vsv2 = vsv2 + step_length*vsv2_kernel
    vsh2 = vsh2 + step_length*vsh2_kernel
    vf2 = vf2 + step_length*vf2_kernel

    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', sqrt(vphi2 + 4.0/3.0*vsh2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', sqrt(vphi2 + 4.0/3.0*vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', sqrt(vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', sqrt(vsh2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', vf2/(vphi2+4.0/3.0*vsh2-2.0*vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho', rho)

  enddo

  !====== Finalize
  if (myrank == 0) close(IOUT)

  call synchronize_all()
  call finalize_mpi()

end program
