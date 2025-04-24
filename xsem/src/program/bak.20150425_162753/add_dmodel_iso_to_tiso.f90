! add_dmodel_iso_to_tiso <mesh_dir> <model_dir> <dmodel_dir> <out_dir> <max_dlnv_allowed> <flag_force_max_dlnv> <log_file>

subroutine selfdoc()
  print *, 'add_dmodel_iso_to_tiso <mesh_dir> <model_dir> <dmodel_dir>&
          & <out_dir> <max_dlnv_allowed> <flag_force_max_dlnv=0/1> <log_file>'
  print *, 'files read: <mesh_dir>/proc000***_reg1_solver_data.bin'
  print *, 'files read: <model_dir>/proc000***_reg1_[vpv,vph,vsv,vsh,eta,rho].bin'
  print *, 'files read: <dmodel_dir>/proc000***_reg1_[mu,lamda_2mu,rho]_dmodel.bin'
  print *, 'files written: <out_dir>/proc000***_reg1_[vpv,vph,vsv,vsh,eta,rho]_new.bin'
  print *, 'max_dlnv_allowed > 0: model + dmodel * step_length'
  print *, 'max_dlnv_allowed < 0: model - dmodel * step_length'
  print *, 'the step_length is determined when, '
  print *, 'flag_force_max_dlnv = 1: force max velocity perturbation = max_dlnv_allowed'
  print *, 'flag_force_max_dlnv = 0: only reduce max velocity perturbation to max_dlnv_allowed'
  print *, '  if already smaller than max_dlnv_allowed, then do nothing'
  print *, 'write(append) to log_file: step_length = '
end subroutine

program sum_kernels_cijkl_to_iso 

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use parallel
  use sem_IO
  use sem_tomography

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 7

  !---- cmd line args
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, dmodel_dir, out_dir, log_file
  double precision :: max_dlnv_allowed
  logical :: flag_force_max_dlnv
  
  !---- local variables
  integer :: sizeprocs, myrank, i, ier
  integer :: ispec
  character(len=MAX_STRING_LEN) :: model_name
  logical :: is_exist ! log_file
  type(sem_mesh) :: mesh_data
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vpv, vph, vsv, vsh, eta, rho
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: A, C, F, L, N
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dmu, dlamda_2mu, drho, dlamda
  double precision :: max_dln_vpv, max_dln_vph, max_dln_vsv, max_dln_vsh, max_dlnv
  double precision :: step_length, max_dlnv_all

  ! ============ program starts here =====================

  ! initializes mpi 
  call initialize()

  !print *, 'myrank=', myrank

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
      call abort_mpi()
    endif
  enddo

  read(args(1),'(a)') mesh_dir 
  read(args(2),'(a)') model_dir 
  read(args(3),'(a)') dmodel_dir 
  read(args(4),'(a)') out_dir 
  read(args(5),*) max_dlnv_allowed
  select case (args(6))
    case ('0') 
      flag_force_max_dlnv = .false.
    case ('1') 
      flag_force_max_dlnv = .true.
    case default
      if (myrank == 0) then
        write(*,*) 'Error: flag_force_max_dlnv must be 0 or 1'
      end if
      call abort_mpi()
  end select
  read(args(7),'(a)') log_file 

  !---- initialize arrays
  call sem_set_dimension(iregion)

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
    dlamda_2mu(NGLLX,NGLLY,NGLLZ,NSPEC), &
           dmu(NGLLX,NGLLY,NGLLZ,NSPEC), &
        dlamda(NGLLX,NGLLY,NGLLZ,NSPEC), &
          drho(NGLLX,NGLLY,NGLLZ,NSPEC))

  ! open log_file for write(append)
  if (myrank == 0) then
    inquire(file=log_file, exist=is_exist)
    if (is_exist) then
      open(unit=IOUT, file=log_file, action='write', &
           position='append', status='old', iostat=ier)
    else
      open(unit=IOUT, file=log_file, action='write', status='new', iostat=ier)
    end if
    if (ier /= 0) then
      print *, 'Error open file for write(append): ', log_file 
      call abort_mpi() 
    end if
    write(IOUT,*)
    write(IOUT,*) '#------------------------------------------#' 
    write(IOUT,*) '# add_dmodel_iso_to_tiso:', trim(timestamp()) 
    write(IOUT,*) '#------------------------------------------#' 
    write(IOUT,*)
  end if

  !read mesh (only use ispec_is_tiso)
  call sem_read_mesh(mesh_data, mesh_dir, myrank, iregion)

  ! read old models
  model_name = 'vpv'; call sem_read_model(vpv, model_dir, myrank, iregion, model_name)
  model_name = 'vph'; call sem_read_model(vph, model_dir, myrank, iregion, model_name)
  model_name = 'vsv'; call sem_read_model(vsv, model_dir, myrank, iregion, model_name)
  model_name = 'vsh'; call sem_read_model(vsh, model_dir, myrank, iregion, model_name)
  model_name = 'eta'; call sem_read_model(eta, model_dir, myrank, iregion, model_name)
  model_name = 'rho'; call sem_read_model(rho, model_dir, myrank, iregion, model_name)
  
  ! enforce isotropy of isotropic elements
  do ispec = 1, NSPEC
    if (mesh_data%ispec_is_tiso(ispec)) then
    else
      vph(:,:,:,ispec) = vpv(:,:,:,ispec)
      vsh(:,:,:,ispec) = vsv(:,:,:,ispec)
      eta(:,:,:,ispec) = 1.0_CUSTOM_REAL
    end if
  end do

! print *,'**** old model:', 'myrank=', myrank
! print *,'min(vpv,vph)=', minval(vpv), minval(vph),' myrank=',myrank
! print *,'max(vpv,vph)=', maxval(vpv), maxval(vph),' myrank=',myrank
! print *,'min(vsv,vsh)=', minval(vsv), minval(vsh),' myrank=',myrank
! print *,'max(vsv,vsh)=', maxval(vsv), maxval(vsh),' myrank=',myrank
! print *,'min,max(rho)=', minval(rho), maxval(rho),' myrank=',myrank
! print *,'min,max(eta)=', minval(eta), maxval(eta),' myrank=',myrank

  ! get A,C,F,L,N from old model
  A = vph**2 * rho; C = vpv**2 * rho;
  L = vsv**2 * rho; N = vsh**2 * rho;
  F = eta * (A - 2*L)

! print *,'min(A,C,F,L,N)=', minval(A),minval(C),minval(F),minval(L),minval(N),' myrank=',myrank
! print *,'max(A,C,F,L,N)=', maxval(A),maxval(C),maxval(F),maxval(L),maxval(N),' myrank=',myrank

  ! read model update step 
  model_name = 'lamda_2mu_dmodel'
  call sem_read_model(dlamda_2mu, dmodel_dir, myrank, iregion, model_name)
  model_name = 'mu_dmodel'
  call sem_read_model(dmu, dmodel_dir, myrank, iregion, model_name)
  model_name = 'rho_dmodel'
  call sem_read_model(drho, dmodel_dir, myrank, iregion, model_name)

! print *,'**** model update direction:', 'myrank=', myrank 
! print *,'min,max(dlamda_2mu)=', minval(dlamda_2mu), maxval(dlamda_2mu),' myrank=',myrank
! print *,'min,max(dmu)=', minval(dmu), maxval(dmu),' myrank=',myrank
! print *,'min,max(drho)=', minval(drho), maxval(drho),' myrank=',myrank

  ! calculate maximum absolute relative perturbation of velocities
  max_dln_vpv = 0.5 * maxval(abs(DBLE(dlamda_2mu)/C - DBLE(drho/rho)))
  max_dln_vph = 0.5 * maxval(abs(DBLE(dlamda_2mu/A - DBLE(drho/rho))))
  max_dln_vsv = 0.5 * maxval(abs(DBLE(dmu)/L - DBLE(drho/rho)))
  max_dln_vsh = 0.5 * maxval(abs(DBLE(dmu)/N - DBLE(drho/rho)))
  max_dlnv = max(max_dln_vpv,max_dln_vph,max_dln_vsv,max_dln_vsh)

  call synchronize_all()
  call max_all_dp(max_dlnv,max_dlnv_all)
  call bcast_all_singledp(max_dlnv_all)

! print *,'**** relative perturbation from one times update direction:', 'myrank=', myrank 
! print *, 'max(abs(dln_vpv))=', max_dln_vpv,' myrank=',myrank
! print *, 'max(abs(dln_vph))=', max_dln_vph,' myrank=',myrank
! print *, 'max(abs(dln_vsv))=', max_dln_vsv,' myrank=',myrank
! print *, 'max(abs(dln_vsh))=', max_dln_vsh,' myrank=',myrank
! print *, 'max(abs(dlnv))=', max_dlnv,' myrank=',myrank

  ! determine the model update step length
  if ((max_dlnv_all>abs(max_dlnv_allowed)) .or. flag_force_max_dlnv) then
    step_length = max_dlnv_allowed / max_dlnv_all
  else ! only set positive/negative
    step_length = SIGN(1.d0, max_dlnv_allowed)
  end if

  if (myrank == 0) then
    write(IOUT,*) 'max_dlnv_all = ', max_dlnv_all
    write(IOUT,*) 'step_length = ', step_length
  end if

! print *, 'step_length=', step_length, ' myrank=',myrank

  ! calculate new models
  dlamda_2mu = SNGL(dlamda_2mu * step_length)
  dmu = SNGL(dmu * step_length)
  drho = SNGL(drho * step_length)
  dlamda = dlamda_2mu - 2*dmu

! print *,'**** model update step:', ' myrank=',myrank
! print *,'min,max(dlamda_2mu)=', minval(dlamda_2mu), maxval(dlamda_2mu),' myrank=',myrank
! print *,'min,max(dmu)=', minval(dmu), maxval(dmu),' myrank=',myrank
! print *,'min,max(dlamda)=', minval(dlamda), maxval(dlamda),' myrank=',myrank
! print *,'min,max(drho)=', minval(drho), maxval(drho),' myrank=',myrank

  rho = rho + drho
  vpv = sqrt((C+dlamda_2mu)/rho)
  vph = sqrt((A+dlamda_2mu)/rho)
  vsv = sqrt((L+dmu)/rho)
  vsh = sqrt((N+dmu)/rho)
  eta = (F+dlamda)/(A-2*L+dlamda)

! print *,'**** new model:', ' myrank=',myrank
! print *,'min(vpv,vph)=', minval(vpv), minval(vph),' myrank=',myrank
! print *,'max(vpv,vph)=', maxval(vpv), maxval(vph),' myrank=',myrank
! print *,'min(vsv,vsh)=', minval(vsv), minval(vsh),' myrank=',myrank
! print *,'max(vsv,vsh)=', maxval(vsv), maxval(vsh),' myrank=',myrank
! print *,'min,max(rho)=', minval(rho), maxval(rho),' myrank=',myrank
! print *,'min,max(eta)=', minval(eta), maxval(eta),' myrank=',myrank

  ! enforce isotropy for isotropic elements
  do ispec = 1, NSPEC
    if (mesh_data%ispec_is_tiso(ispec)) then
      ! todo: why if (.not. xxx) does not work?
    else
      vph(:,:,:,ispec) = vpv(:,:,:,ispec)
      vsh(:,:,:,ispec) = vsv(:,:,:,ispec)
      eta(:,:,:,ispec) = 1.0_CUSTOM_REAL
    end if
  end do

  ! write new models
  model_name = 'vpv_new'; call sem_write_model(vpv, out_dir, myrank, iregion, model_name)
  model_name = 'vph_new'; call sem_write_model(vph, out_dir, myrank, iregion, model_name)
  model_name = 'vsv_new'; call sem_write_model(vsv, out_dir, myrank, iregion, model_name)
  model_name = 'vsh_new'; call sem_write_model(vsh, out_dir, myrank, iregion, model_name)
  model_name = 'eta_new'; call sem_write_model(eta, out_dir, myrank, iregion, model_name)
  model_name = 'rho_new'; call sem_write_model(rho, out_dir, myrank, iregion, model_name)

  ! stop all the MPI processes, and exit
  if (myrank == 0) close(IOUT)
  call synchronize_all()
  call finalize_mpi()

!=============================================================================
contains

subroutine initialize()

  implicit none

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (sizeprocs /= NPROCTOT_VAL) then
    write(*,*) 'Error detected, aborting MPI... proc ',myrank
    if (myrank == 0) then
      print*, 'Error number of processors supposed to run on : ',NPROCTOT_VAL
      print*, 'Error number of MPI processors actually run on: ',sizeprocs
      print*
      print*, 'please rerun with: mpirun -np ',NPROCTOT_VAL,' bin/xadd_model .. '
      call selfdoc()
    end if
    call abort_mpi()
  endif

end subroutine initialize

end program sum_kernels_cijkl_to_iso
