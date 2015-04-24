!
! interpolate model values from a list of (x,y,z) values
!

program interp_xyz

  use constants, only: NPROCTOT_VAL
  use mesh

  implicit none

  character(len=1), parameter :: delimiter = ','
  integer, parameter :: MAX_NUM_XYZ = 100000
  integer, parameter :: MAX_NUM_MODEL = 10

  ! local variables
  character(len=MAX_STRING_LEN) :: arg(5)

  integer :: i,j, ier, iproc, num_xyz
  character(len=MAX_STRING_LEN) :: xyz_list, model_str, out_file
  character(len=255) :: strtok

  real(RK) :: xyz_interp(3,MAX_NUM_XYZ)
  real(RK), allocatable :: res_dist(:), model_interp(:,:)
  integer, allocatable :: iloc(:)

  ! read command line arguments
  do i = 1, 5
    call get_command_argument(i,arg(i), status=ier)
    if (i <= 5 .and. trim(arg(i)) == '') then
      print *, 'Usage: interp_xyz xyz_list model_str model_dir topo_dir out_file'
      stop
    endif
  enddo
  read(arg(1),'(a)') xyz_list
  read(arg(2),'(a)') model_str
  read(arg(3),'(a)') model_dir
  read(arg(4),'(a)') topo_dir
  read(arg(5),'(a)') out_file
  iregion = 1

  ! read xyz list
  open(unit=IIN, file=xyz_list, status='old', iostat=ier)
  if (ier /= 0) then
    print *, 'Error open file: ', xyz_list 
    stop
  endif
  num_xyz = 0
  do
    read(IIN,*, iostat=ier) xyz_interp(:,num_xyz+1)
    if (ier /= 0) exit
    num_xyz = num_xyz + 1
    if (num_xyz >= MAX_NUM_XYZ) exit
  enddo
  close(IIN)

  print *, 'num_xyz=', num_xyz
 
  ! parse model names
  allocate(model_names(MAX_NUM_MODEL))
  i = 1
  model_names(i) = trim(strtok(model_str,delimiter))
  do while (model_names(i) /= char(0))
    i = i+1
    if (i > MAX_NUM_MODEL) exit
    model_names(i) = trim(strtok(char(0),delimiter))
  enddo
  num_model = i-1

  print *, 'num_model=', num_model

  ! do interpolation in each mesh slice
  allocate(iloc(num_xyz), res_dist(num_xyz), &
           model_interp(num_xyz,num_model))
  res_dist = HUGEVAL_SNGL
  iloc = -1
  model_interp(:,:) = 0.0_RK

  do iproc = 1, NPROCTOT_VAL
  !do iproc = 1, 1
    !print *, 'slice ', iproc
    myrank = iproc - 1
    call read_mesh()
    call read_model()
    call interp_model(xyz_interp(:,1:num_xyz), model_interp, iloc, res_dist)
  enddo

  ! write out results
  open(unit=IOUT, file=out_file, status='unknown', iostat=ier)
  if (ier /= 0) then
    print *, 'Error open file: ', out_file
    stop
  endif
  do i = 1, num_xyz
    write(IOUT,'($,3E15.7,a,I2,2x,E15.7,a)') xyz_interp(:,i),' | ', &
      iloc(i),res_dist(i)/SNGL(MESH_SIZE),' | '
    do j = 1, num_model
      write(IOUT,'($,E15.7)') model_interp(i,j)
    enddo
    write(IOUT,'()')
  enddo
  close(IOUT)

end program
