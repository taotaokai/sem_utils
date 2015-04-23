! sem_inner_product <mesh_dir> <model_dir_1> <model_name_1> 
!                              <model_dir_2> <model_name_2> 
!                              <log_file>    <key_name>

subroutine selfdoc()
  print *, 'Usage: sem_inner_product <mesh_dir> <model_dir_1> <model_name_1> &
    & <model_dir_2> <model_name_2> <log_file> <key_name>'
  print *, 'compute the volume integral of (model_1 * model_2)'
  print *, 'model_name: vpv, rho_model, mu_model etc.'
  print *, 'model file read: <model_dir>/proc000***_reg1_<model_name>.bin'
  print *, 'append to log_file: <key_name> = <volume integral result>'
  stop
end subroutine

program sem_inner_product

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use sem_IO
  use sem_tomography

  implicit none

  !---- parameters 
  integer, parameter :: iregion = 1
  integer, parameter :: nargs = 7

  !---- cmd line args
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir_1, model_name_1 
  character(len=MAX_STRING_LEN) :: model_dir_2, model_name_2 
  character(len=MAX_STRING_LEN) :: log_file, key_name

  !---- local variables
  integer :: iproc, i, ier
  logical :: is_exist ! log_file
  type(mesh) :: mesh_data
  real(CUSTOM_REAL), dimension(:,:,:,:), allocatable :: model_1, model_2
  real(kind=8), dimension(:,:,:,:), allocatable :: wgll_vol, model_prod
  real(kind=8) :: summation, volume_integral

  ! ============ program starts here =====================

  !---- get command line arguments 
  do i = 1, nargs
    call get_command_argument(i,args(i), status=ier)
    if (trim(args(i)) == '') then
      call selfdoc()
    endif
  enddo

  read(args(1),'(a)') mesh_dir
  read(args(2),'(a)') model_dir_1
  read(args(3),'(a)') model_name_1
  read(args(4),'(a)') model_dir_2
  read(args(5),'(a)') model_name_2
  read(args(6),'(a)') log_file
  read(args(7),'(a)') key_name

  !---- initialize arrays
  call sem_set_dimension(iregion)
  allocate( model_1(NGLLX,NGLLY,NGLLZ,NSPEC), &
            model_2(NGLLX,NGLLY,NGLLZ,NSPEC), &
           wgll_vol(NGLLX,NGLLY,NGLLZ,NSPEC), &
         model_prod(NGLLX,NGLLY,NGLLZ,NSPEC))

  call sem_init_mesh(mesh_data)
  
  !---- operate on each slice 
  volume_integral = 0.d0
  do iproc = 0, NPROCTOT_VAL-1

    ! read mesh
    call sem_read_mesh(mesh_data, mesh_dir, iproc, iregion)

    ! read models
    call sem_read_model(model_1, model_dir_1, iproc, iregion, model_name_1)
    call sem_read_model(model_2, model_dir_2, iproc, iregion, model_name_2)

    ! inner product(model_1,model_2)
    call sem_volume_weight_gll(mesh_data, wgll_vol)

    model_prod = wgll_vol
    model_prod = model_prod * model_1
    model_prod = model_prod * model_2
    call Kahan_sum(model_prod, summation)

    ! accumulate result
    volume_integral = volume_integral + summation

  end do

  !---- write out log file
  inquire(file=log_file, exist=is_exist)
  if (is_exist) then
    open(unit=IOUT, file=log_file, action='write', &
         position='append', status='old', iostat=ier)
  else
    open(unit=IOUT, file=log_file, action='write', status='new', iostat=ier)
  end if
  if (ier /= 0) then
    print *, 'Error open file for write(append): ', log_file 
    stop
  end if
  write(IOUT,*)
  write(IOUT,*) '#------------------------------------------#' 
  write(IOUT,'("# (",a,",",a,"): ",a)') &
    trim(model_name_1), trim(model_name_2), trim(timestamp()) 
  write(IOUT,*) '#------------------------------------------#' 
  write(IOUT,*)
  write(IOUT,'(a," = ",SP,ES15.7)') trim(key_name), volume_integral 
  close(IOUT)

! =========================
contains

subroutine Kahan_sum(model_array, model_sum)
  
  real(kind=8), dimension(:,:,:,:), intent(in) :: model_array
  real(kind=8), intent(out) :: model_sum

  integer :: i,j,k,ispec
  real(kind=8) :: t, y, c

  c = 0.d0
  model_sum = 0.d0

  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          y = model_array(i,j,k,ispec) - c
          t = model_sum + y
          c = (t - model_sum) - y
          model_sum = t
        end do
      end do
    end do
  end do

end subroutine

end program
