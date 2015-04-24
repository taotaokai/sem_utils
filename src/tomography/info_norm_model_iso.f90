subroutine selfdoc()
  write(*,'(a)') &

'NAME', &
'', &
'   info_model_norm - print norm of model files', &
'', &
'SYNOPSIS', &
'', &
'  info_model_norm <mesh_dir> <model_dir> <model_names>', &
'', &
'DESCRIPTION', &
'', &
'  norm = sqrt(sum(volume integrals of each model_name))', &
'', &
'COMMENT', &
'', &
'  The mesh topology should be the same for all models.', &
'', &
'PARAMETERS', &
'', &
'  <mesh_dir>  proc000***_reg1_solver_data.bin', &
'  <model_dir>  proc000***_reg1_{model_name}.bin', &
'  <model_names> comma seperated model names, e.g. mu_kernel,rho_kernel', &
''
end subroutine

program info_model_norm

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use sem_IO
  use sem_tomography

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 3

  !---- cmd line args
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir
  
  !---- local variables
  integer :: i, ier
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: model_names
  character(len=1024) :: dummy
  type(sem_mesh) :: mesh_data
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: model_array
  real(kind=CUSTOM_REAL) :: norm2_sum

  ! ============ program starts here =====================
  
  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
      stop
    endif
  enddo

  read(args(1),'(a)') mesh_dir 
  read(args(2),'(a)') model_dir
  read(args(3),'(a)') model_names 

  allocate(dmodel_dir(nstep_lbfgs), dkernel_dir(nstep_lbfgs))
  allocate(step_length(nstep_lbfgs))

  rewind(IIN)
  do i = 1, nstep_lbfgs
    read(IIN,*,iostat=ier) dmodel_dir(i), dkernel_dir(i), step_length(i)
    if (ier /= 0) then
       write(*,*) 'Error read line number ', i
       stop
    endif
  enddo
  close(IIN)

  !---- initialize data 
  call sem_set_dimension(iregion)

  allocate(rho(nstep_lbfgs), alpha(nstep_lbfgs))
  allocate( dmodel_old(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC,nstep_lbfgs), &
           dkernel_old(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC,nstep_lbfgs) )
  allocate(dmodel_new(nmodel,NGLLX,NGLLY,NGLLZ,NSPEC))

  ! get mesh data
  call sem_read_mesh(mesh_data, mesh_dir, myrank, iregion)
  call sem_set_volume_gll(mesh_data)

  ! set model names
  allocate(kernel_names(3), dmodel_names(3), dkernel_names(3))

  kernel_names(1)  = 'lamda_2mu_kernel'
  dmodel_names(1)  = 'lamda_2mu_dmodel'
  dkernel_names(1) = 'lamda_2mu_dkernel'

  kernel_names(2)  = 'mu_kernel'
  dmodel_names(2)  = 'mu_dmodel'
  dkernel_names(2) = 'mu_dkernel'

  kernel_names(3)  = 'rho_kernel'
  dmodel_names(3)  = 'rho_dmodel'
  dkernel_names(3) = 'rho_dkernel'

  !---- two-loop l-bfgs

  ! read kernel in the current iteration: g0
  call sem_read_model_n(dmodel_new, kernel_dir, myrank, iregion, nmodel, kernel_names)

  ! q = - g0
  dmodel_new = -1.0_CUSTOM_REAL * dmodel_new

  ! first loop: i -> i'th iteration before the current
  if (myrank == 0) then
    write(IOUT,*) '#---- l-FBGS: first loop'
  end if

  do i = 1, nstep_lbfgs

    if (myrank==0) then
      write(IOUT,*) '# i=',i, ' -th iteration before'
    end if

    ! read old dm(i), dg(i)
    call sem_read_model_n(dmodel_old(:,:,:,:,:,i), dmodel_dir(i), &
      myrank, iregion, nmodel, dmodel_names)
    call sem_read_model_n(dkernel_old(:,:,:,:,:,i), dkernel_dir(i), &
      myrank, iregion, nmodel, dkernel_names)

    ! dm(i) = step_length * dm(i)
    dmodel_old(:,:,:,:,:,i) = step_length(i) * dmodel_old(:,:,:,:,:,i)

    ! rho(i) = 1 / (dm(i), dg(i))
    call sem_inner_prod(mesh_data, nmodel, dmodel_old(:,:,:,:,:,i), &
                        dkernel_old(:,:,:,:,:,i), rho(i))
    call synchronize_all()
    call sum_all_cr(rho(i), rho_sum)
    call bcast_all_singlecr(rho_sum)
    rho(i) = 1.0_CUSTOM_REAL / rho_sum

    if (myrank == 0) then
       write(IOUT,*) 'rho(i)=1/(dm(i),dg(i))=', rho(i)
    end if

    ! alpha(i) = rho(i) * (dm(i), q)
    call sem_inner_prod(mesh_data, nmodel, dmodel_old(:,:,:,:,:,i), &
                        dmodel_new, alpha(i))
    ! gather all alpha(i)
    call synchronize_all()
    call sum_all_cr(alpha(i), alpha_sum)
    call bcast_all_singlecr(alpha_sum)
    alpha(i) = rho(i) * alpha_sum

    if (myrank == 0) then
      write(IOUT,*) 'alpha(i)=rho(i)*(dm(i),q)=', alpha(i)
    end if

    ! q = q - alpha(i)*dg(i)
    dmodel_new = dmodel_new - alpha(i)*dkernel_old(:,:,:,:,:,i)
  end do

  ! initial Hessian: q = H0 * q
  ! H0 ~ (dm(1), dg(1)) / (dg(1), dg(1))
  call sem_inner_prod(mesh_data, nmodel, dkernel_old(:,:,:,:,:,1), &
                      dkernel_old(:,:,:,:,:,1), gamma)
  call synchronize_all()
  call sum_all_cr(gamma, gamma_sum)
  call bcast_all_singlecr(gamma_sum)
  gamma = 1.0_CUSTOM_REAL / rho(1) / gamma_sum
  dmodel_new = gamma * dmodel_new

  if (myrank == 0) then
    write(IOUT,*) '#---- apply initial Hessian'
    write(IOUT,*) '# H0 ~ gamma * I ~ (dm(1), dg(1)) / (dg(1), dg(1)) * I'
    write(IOUT,*) 'gamma=', gamma
  end if

  ! second loop of lbfgs
  if (myrank==0) then 
    write(IOUT,*) '#---- l-FBGS, second loop ****'
  end if

  do i = nstep_lbfgs, 1, -1

    if (myrank==0) then
      write(IOUT,*) '# i=',i, ' -th iteration before'
    end if

    ! beta = rho(i) * (dg(i), q)
    call sem_inner_prod(mesh_data, nmodel, dmodel_new, dkernel_old(:,:,:,:,:,i), beta)
    call synchronize_all()
    call sum_all_cr(beta, beta_sum)
    call bcast_all_singlecr(beta_sum)
    beta = rho(i) * beta_sum

    if (myrank == 0) then 
      write(IOUT,*) 'beta=rho(i)*(dg(i),q)=', beta
    end if

    ! q = q + (alpha(i) - beta)*dm(i)
    dmodel_new = dmodel_new + (alpha(i) - beta) * dmodel_old(:,:,:,:,:,i)
  end do

  !---- write out new dmodel
  call synchronize_all()
  call sem_write_model_n(dmodel_new, out_dir, myrank, iregion, nmodel, dmodel_names)

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
      print*, 'please rerun with: mpirun -np ',NPROCTOT_VAL
      print*
      call selfdoc()
    end if
    call abort_mpi()
  endif

end subroutine initialize

end program