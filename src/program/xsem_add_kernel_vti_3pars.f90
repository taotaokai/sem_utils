subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_kernel_vti_3pars"
  print '(a)', "    - add vti (precondioned-)kernel (vp2, vsv2, vsh2)"
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
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_[vp2,vsv2,vsh2]_kernel.bin"
  print '(a)', "  (string) kernel_suffix:  proc*_reg1_[vp2,vsv2,vsh2]_kernel<kernel_suffix>.bin"
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
  print '(a)', "  4. thomsen's anisotropic parameters: eps = delta = 0 => C11 = C33, C13 = C11 - 2*C44"
  print '(a)', "     vp2_kernel -> C11,C33, vsv2_kernel -> C44, vsh2_kernel -> C66"
  print '(a)', "  5. For isotropic element the program will reduce the vti kernel to isotropic kernel and alos enforce the output be isotropic."

end subroutine


program xsem_add_kernel_vti_3pars

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
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vph2, vpv2, vsv2, vsh2, vf2, rho
  ! kernel
  real(dp), dimension(:,:,:,:), allocatable :: vp2_kernel, vsv2_kernel, vsh2_kernel
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
        print *, '[ERROR] flag_force_max_dlnv must be 0 or 1'
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
      write(*,*) '[ERROR] failed to open file ', trim(log_file)
      call abort_mpi()
    endif

    write(IOUT, '(a)') '#[LOG] xsem_add_kernel_vti_3pars'

  endif

  !===== Get step length first

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(vph2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vf2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vp2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dln_rho(NGLLX,NGLLY,NGLLZ,nspec))
 
  ! get maximum velocity perturbation from dmodel
  max_dlnv = 0.0_dp

  if (myrank ==0) write(IOUT,'(a)') "#====== get step_length"
  call synchronize_all()

  do iproc = myrank, (nproc-1), nrank

    ! read in mesh 
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)

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

    ! check if element is isotropic
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        ! enforce model to be isotropic
        vsh2(:,:,:,ispec) = vsv2(:,:,:,ispec)
        ! reduce vti kernel to iso kernel
        vsv2_kernel(:,:,:,ispec) = vsv2_kernel(:,:,:,ispec) + vsh2_kernel(:,:,:,ispec)
        vsh2_kernel(:,:,:,ispec) = vsv2_kernel(:,:,:,ispec)
      endif
    enddo

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

    ! read in mesh 
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', vf2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    ! convert to squared value and also vf2
    vph2 = vph2**2              ! = C11/rho
    vpv2 = vpv2**2              ! = C33/rho
    vsv2 = vsv2**2              ! = C44/rho
    vsh2 = vsh2**2              ! = C66/rho
    vf2 = vf2*(vph2 - 2.0*vsv2) ! = C13/rho

    ! read kernels
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vp2_kernel'//kernel_suffix, vp2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsv2_kernel'//kernel_suffix, vsv2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsh2_kernel'//kernel_suffix, vsh2_kernel)

    ! check if element is isotropic
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        ! enforce model to be isotropic
        vph2(:,:,:,ispec) = vpv2(:,:,:,ispec)
        vsh2(:,:,:,ispec) = vsv2(:,:,:,ispec)
        vf2(:,:,:,ispec) = vpv2(:,:,:,ispec) - 2.0*vsv2(:,:,:,ispec)
        ! reduce vti kernel to iso kernel
        vsv2_kernel(:,:,:,ispec) = vsv2_kernel(:,:,:,ispec) + vsh2_kernel(:,:,:,ispec)
        vsh2_kernel(:,:,:,ispec) = vsv2_kernel(:,:,:,ispec)
      endif
    enddo

    ! scale
    ! voigt_G = (2*C44 + C66)/3
    dln_rho = ratio_rho*step_length*(2.0*vsv2_kernel + vsh2_kernel)/(2.0*vsv2 + vsh2)/2.0
    rho = rho*(1.0 + dln_rho)
    print *, "dln_rho: min/max = ", minval(dln_rho), maxval(dln_rho)

    vph2 = vph2 + step_length*vp2_kernel
    vpv2 = vpv2 + step_length*vp2_kernel
    vsv2 = vsv2 + step_length*vsv2_kernel
    vsh2 = vsh2 + step_length*vsh2_kernel
    vf2  = vf2  + step_length*(vp2_kernel - 2.0*vsv2_kernel)

    ! write new models
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', sqrt(vph2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', sqrt(vpv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', sqrt(vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', sqrt(vsh2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', vf2/(vph2 - 2.0*vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho', rho)

  enddo

  !====== Finalize
  if (myrank == 0) close(IOUT)

  call synchronize_all()
  call finalize_mpi()

end program
