subroutine selfdoc()

  write(*,'(a)') &

'NAME', &
'  get_dmodel_iso_lbfgs <mesh_dir> <dmodel_dkernel_list> <kernel_dir> <out_dir>', &
'', &
'DESCRIPTION', &
'  get model update from the kernel of this iteration and the (dmodel,dkernel)', &
'from the previous iterations using l-BFGS method', &
'', &
'SYNOPSIS', &
'  get_dmodel_iso_lbfgs <mesh_dir> <dmodel_dkernel_list> <out_dir>', &
'', &
'COMMENT', &
'  The mesh geometry for all the iterations should be the same.', &
'', &
'PARAMETERS', &
'  <mesh_dir>  proc000***_reg1_solver_data.bin', &
'  <dmodel_dkernel_list>  a text file containing lines of ', &
'    dMODEL_dir dKERNEL_dir step_length (dm,dg) (dg,dg) (dm,dm)', &
'  <kernel_dir>  proc000***_reg1_[mu,lamda_2mu,rho]_kernel.bin', &
'  <out_dir>  proc000***_reg1_[mu,lamda_2mu,rho]_dmodel.bin'

end subroutine

program get_dmodel_iso_lbfgs

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use parallel
  use sem_IO
  use sem_tomography

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 4

  !---- cmd line args
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, dmodel_dkernel_list, kernel_dir, out_dir
  
  !---- local variables
  integer :: i, ier, sizeprocs, myrank
  integer :: ispec, nstep_lbfgs
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: dmodel_dir, dkernel_dir
  character(len=1024) :: dummy
  type(mesh) :: mesh_data
  type(model_iso) :: dmodel_iso_new
  type(model_iso), dimension(:), allocatable :: dmodel_iso, dkernel_iso
  double precision, dimension(:), allocatable :: step_length, rho, alpha
  double precision :: beta, rho_sum, alpha_sum, beta_sum, gamma, gamma_sum
  real(kind=8), dimension(:,:,:,:), allocatable :: wgll_vol

  ! ============ program starts here =====================
  
  ! initializes mpi 
  call initialize()

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
      call abort_mpi()
    endif
  enddo

  read(args(1),'(a)') mesh_dir 
  read(args(2),'(a)') dmodel_dkernel_list
  read(args(3),'(a)') kernel_dir
  read(args(4),'(a)') out_dir 

  !---- read dmodel_dkernel_list
  open(unit=IIN, file=trim(dmodel_dkernel_list), status='old', action='read', iostat=ier)
  if (ier /= 0) then
     write(*,*) 'Error open file: ', trim(dmodel_dkernel_list)
     stop
  endif
  nstep_lbfgs = 0 ! get number of lines
  do
    read(IIN,'(a)',iostat=ier) dummy
    if (ier /= 0) exit
    nstep_lbfgs = nstep_lbfgs + 1
  enddo

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

  !---- initialize arrays
  call sem_set_dimension(iregion)

  allocate(rho(nstep_lbfgs), alpha(nstep_lbfgs), &
           wgll_vol(NGLLX,NGLLY,NGLLZ,NSPEC))
  allocate(dmodel_iso(nstep_lbfgs), dkernel_iso(nstep_lbfgs))

  call sem_init_mesh(mesh_data)
 
  !---- get volume integration weight
  call sem_read_mesh(mesh_data, mesh_dir, myrank, iregion)
  call sem_volume_weight_gll(mesh_data, wgll_vol)

  !---- two-loop l-bfgs 
  
  ! read kernel
  call sem_read_model_iso(dmodel_iso_new, kernel_dir, myrank, iregion)

  ! q = - g
  dmodel_iso_new%lamda_2mu = -1.0_CUSTOM_REAL * dmodel_iso_new%lamda_2mu
  dmodel_iso_new%mu        = -1.0_CUSTOM_REAL * dmodel_iso_new%mu
  dmodel_iso_new%rho       = -1.0_CUSTOM_REAL * dmodel_iso_new%rho

  ! first loop
  do i = 1, nstep_lbfgs

    ! read dm(i), dg(i)
    call sem_read_model_iso(dmodel_iso(i), dmodel_dir(i), myrank, iregion)
    call sem_read_model_iso(dkernel_iso(i), dkernel_dir(i), myrank, iregion)
    ! dm(i) = step_length * dm(i)
    dmodel_iso(i)%lamda_2mu = step_length(i) * dmodel_iso(i)%lamda_2mu
    dmodel_iso(i)%mu        = step_length(i) * dmodel_iso(i)%mu
    dmodel_iso(i)%rho       = step_length(i) * dmodel_iso(i)%rho

    ! rho(i) = 1 / (dm(i), dg(i))
    call sem_inner_product_iso(dmodel_iso(i), dkernel_iso(i), wgll_vol, rho(i))
    call synchronize_all()
    call sum_all_dp(rho(i), rho_sum)
    call bcast_all_singledp(rho_sum)
    rho(i) = 1.0d0 / rho_sum

    ! alpha(i) = rho(i) * (dm(i), q)
    call sem_inner_product_iso(dmodel_iso(i), model_iso, wgll_vol, alpha(i))
    ! gather all alpha(i)
    call synchronize_all()
    call sum_all_dp(alpha(i), alpha_sum)
    call bcast_all_singledp(alpha_sum)
    alpha(i) = rho(i) * alpha_sum

    ! q = q - alpha(i)*dg(i)
    dmodel_iso_new%lamda_2mu = dmodel_iso_new%lamda_2mu - alpha(i)*dkernel_iso(i)%lamda_2mu
    dmodel_iso_new%mu        = dmodel_iso_new%mu        - alpha(i)*dkernel_iso(i)%mu
    dmodel_iso_new%rho       = dmodel_iso_new%rho       - alpha(i)*dkernel_iso(i)%rho
  end do

  ! initial Hessian: q = H0 * q
  ! H0 ~ (dm(1), dg(1)) / (dg(1), dg(1))
  call sem_inner_product_iso(dkernel_iso(1), dkernel_iso(1), wgll_vol, gamma)
  call synchronize_all()
  call sum_all_dp(gamma, gamma_sum)
  call bcast_all_singledp(gamma_sum)
  gamma = 1.0d0 / rho(1) / gamma_sum
  dmodel_iso_new%lamda_2mu = gamma * dmodel_iso_new%lamda_2mu 
  dmodel_iso_new%mu        = gamma * dmodel_iso_new%mu 
  dmodel_iso_new%rho       = gamma * dmodel_iso_new%rho

  ! second loop
  do i = nstep_lbfgs, 1, -1

    ! beta = rho(i) * (dg(i), q)
    call sem_inner_product_iso(model_iso, dkernel_iso(i), wgll_vol, beta)
    call synchronize_all()
    call sum_all_dp(beta, beta_sum)
    call bcast_all_singledp(beta_sum)
    beta = rho(i) * beta_sum

    ! q = q + (alpha(i) - beta)*dm(i)
    dmodel_iso_new%lamda_2mu  = dmodel_iso_new%lamda_2mu  + (alpha(i) - beta) * dmodel_iso(i)%lamda_2mu 
    dmodel_iso_new%mu         = dmodel_iso_new%mu         + (alpha(i) - beta) * dmodel_iso(i)%mu        
    dmodel_iso_new%rho        = dmodel_iso_new%rho        + (alpha(i) - beta) * dmodel_iso(i)%rho       
  end do

  !---- write out new dmodel
  call synchronize_all()
  call sem_write_model_iso(model_iso, out_dir, myrank, iregion)

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

end program