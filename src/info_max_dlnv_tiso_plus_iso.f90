subroutine selfdoc()
!  write(*,'(a)') &
  print *, 'info_max_dlnv_tiso_plus_iso <model_tiso_dir> <dmodel_iso_dir>'
  print *, 'files read: <model_dir>/proc000***_reg1_[vpv,vph,vsv,vsh,eta,rho].bin'
  print *, 'files read: <dmodel_dir>/proc000***_reg1_[mu,lamda_2mu,rho]_dmodel.bin'
end subroutine

program info_max_dlnv_tiso_plus_iso

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use sem_IO
  use sem_tomography

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 2

  !---- cmd line args
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: model_dir, dmodel_dir
  
  !---- local variables
  integer :: myrank,i, ier
  character(len=MAX_STRING_LEN) :: model_name
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vpv, vph, vsv, vsh, eta, rho
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: A, C, F, L, N
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dmu, dlamda_2mu, drho, dlamda
  double precision, dimension(NPROCTOT_VAL) :: min_dln_vpv, min_dln_vph, min_dln_vsv, min_dln_vsh
  double precision, dimension(NPROCTOT_VAL) :: max_dln_vpv, max_dln_vph, max_dln_vsv, max_dln_vsh

  ! ============ program starts here =====================

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
      stop
    endif
  enddo

  read(args(1),'(a)') model_dir 
  read(args(2),'(a)') dmodel_dir 

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

  !---- loop each slice
  do i = 1, NPROCTOT_VAL

    myrank = i-1

    ! read old models
    model_name = 'vpv'; call sem_read_model(vpv, model_dir, myrank, iregion, model_name)
    model_name = 'vph'; call sem_read_model(vph, model_dir, myrank, iregion, model_name)
    model_name = 'vsv'; call sem_read_model(vsv, model_dir, myrank, iregion, model_name)
    model_name = 'vsh'; call sem_read_model(vsh, model_dir, myrank, iregion, model_name)
    model_name = 'eta'; call sem_read_model(eta, model_dir, myrank, iregion, model_name)
    model_name = 'rho'; call sem_read_model(rho, model_dir, myrank, iregion, model_name)
    
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
  
    ! read dmodel (model update)
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
    min_dln_vpv(i) = 0.5 * minval(DBLE(dlamda_2mu)/C - DBLE(drho/rho))
    max_dln_vpv(i) = 0.5 * maxval(DBLE(dlamda_2mu)/C - DBLE(drho/rho))

    min_dln_vph(i) = 0.5 * minval(DBLE(dlamda_2mu/A - DBLE(drho/rho)))
    max_dln_vph(i) = 0.5 * maxval(DBLE(dlamda_2mu/A - DBLE(drho/rho)))

    min_dln_vsv(i) = 0.5 * minval(DBLE(dmu)/L - DBLE(drho/rho))
    max_dln_vsv(i) = 0.5 * maxval(DBLE(dmu)/L - DBLE(drho/rho))

    min_dln_vsh(i) = 0.5 * minval(DBLE(dmu)/N - DBLE(drho/rho))
    max_dln_vsh(i) = 0.5 * maxval(DBLE(dmu)/N - DBLE(drho/rho))

    print *, '==== myrank=', myrank
    print *, '  dln_vpv: min/max=', min_dln_vpv(i), max_dln_vpv(i)
    print *, '  dln_vph: min/max=', min_dln_vph(i), max_dln_vph(i) 
    print *, '  dln_vsv: min/max=', min_dln_vsv(i), max_dln_vsv(i)
    print *, '  dln_vsh: min/max=', min_dln_vsh(i), max_dln_vsh(i)

  end do ! myrank
  
  print *, '==== statistics'
  print *, '  dln_vpv: min/max=', minval(min_dln_vpv), maxval(max_dln_vpv)
  print *, '  dln_vph: min/max=', minval(min_dln_vph), maxval(max_dln_vph) 
  print *, '  dln_vsv: min/max=', minval(min_dln_vsv), maxval(max_dln_vsv)
  print *, '  dln_vsh: min/max=', minval(min_dln_vsh), maxval(max_dln_vsh)

end program info_max_dlnv_tiso_plus_iso
