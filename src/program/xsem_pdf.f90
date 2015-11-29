subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_pdf "
  print '(a)', "    -  compute PDF of volumetric integral of model amplitude"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_statis \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_name> <nbin> "
  print '(a)', "    <bin_amplitude> <out_file>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_name>.bin"
  print '(a)', "  (string) model_name:  model name, e.g. mu_kernel "
  print '(a)', "  (int) nbin:  number of bins"
  print '(a)', "  (int) use_abs:  flag to bin absolute value |z|, must be 0 or 1"
  print '(a)', "  (string) out_file:  output file name "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. output is a list of <z0> <z1> <percentage_of_volumetric_amplitude_integral>"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_pdf
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir, model_name
  integer :: nbin
  logical :: use_abs
  character(len=MAX_STRING_LEN) :: out_file

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  real(dp), allocatable :: gll_volume(:,:,:,:)
  ! model
  real(dp), allocatable :: model(:,:,:,:)
  ! PDF
  real(dp) :: zmin, zmax, min_i, max_i, bin_size
  real(dp) :: zmin_all, zmax_all
  real(dp), allocatable :: z(:), pdf(:), pdf_all(:)
  real(dp) :: pdf_all_sum, volume, volume_sum

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_pdf: check your input arguments."
      call abort_mpi()
    endif 
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_name
  read(args(5), *) nbin
  select case (args(6))
    case ('0')
      use_abs = .false.
    case ('1')
      use_abs = .true.
    case default
      if (myrank==0) then
        print *, '[ERROR]: use_abs must be 0 or 1'
        call abort_mpi()
      endif
  end select
  read(args(7), '(a)') out_file

  !===== loop each mesh/model slice

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize arrays
  allocate(model(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll_volume(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(z(0:nbin), pdf(nbin))
  if (myrank == 0) allocate(pdf_all(nbin))

  ! get model value range
  zmax = -1.0 * huge(1.0_dp)
  zmin = huge(1.0_dp)
  do iproc = myrank, (nproc-1), nrank
    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    max_i = maxval(model)
    min_i = minval(model)
    if (max_i > zmax) zmax = max_i
    if (min_i < zmin) zmin = min_i
  enddo ! iproc

  call synchronize_all()
   
  call max_all_dp(zmax, zmax_all)
  call min_all_dp(zmin, zmin_all)
  call bcast_all_singledp(zmax_all)
  call bcast_all_singledp(zmin_all)

  call synchronize_all()

  if (myrank == 0) print *, '# zmin/zmax=', zmin_all, zmax_all

  ! compute PDF
  zmax = max(abs(zmax_all), abs(zmin_all))
  if (use_abs) then
    ! bin model amplitude: (0, max|z|)
    bin_size = zmax / nbin
    z = (/(i*bin_size, i=0,nbin)/)
  else
    ! symmetric about 0: (-max|z|, max|z|)
    bin_size = 2.0 * zmax / nbin
    z = (/(-zmax + i*bin_size, i=0,nbin)/)
  endif
  pdf = 0.0_dp
  volume = 0.0_dp
  do iproc = myrank, (nproc-1), nrank
    ! read model
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model)
    ! get gll_volume
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_mesh_gll_volume(mesh_data, gll_volume)
    ! get volumetric integral of model amplitudes within each bin 
    if (use_abs) model = abs(model)
    do i = 1, nbin
      pdf(i) = pdf(i) + sum(gll_volume*abs(model), &
        mask=(model>=z(i-1) .and. model<z(i)) )
    enddo
    ! get volume
    volume = volume + sum(gll_volume)
  enddo
  
  call synchronize_all()

  call sum_all_dp_n(pdf, pdf_all, nbin)
  call sum_all_dp(volume, volume_sum)

  call synchronize_all()

  ! normalize PDF
  if (myrank == 0) then
    pdf_all_sum = sum(pdf_all)
    pdf_all = pdf_all / pdf_all_sum
  endif

  ! write out results
  if (myrank == 0) then
    open(IOUT, file=trim(out_file), status='unknown', &
      form='formatted', action='write', iostat=ier)

    if (ier /= 0) then
      write(*,*) '[ERROR] xsem_pdf: failed to open file ', trim(out_file)
      call abort_mpi() 
    endif

    write(IOUT, '(a)') "# PDF of gll model"
    write(IOUT, '(a,2X,I4)') "# nproc= ", nproc
    write(IOUT, '(a,2X,a)') "# mesh_dir= ", trim(mesh_dir)
    write(IOUT, '(a,2X,a)') "# model_dir= ", trim(model_dir)
    write(IOUT, '(a,2X,a)') "# model_name= ", trim(model_name)
    write(IOUT, '(a,2X,I4)') "# nbin= ", nbin
    write(IOUT, '(a,2X,L2)') "# use_abs= ", use_abs
    write(IOUT, '(a,2X,E12.4)') "# min(z)= ", zmin_all 
    write(IOUT, '(a,2X,E12.4)') "# max(z)= ", zmax_all 
    write(IOUT, '(a,2X,E12.4)') "# volume= ", volume 
    write(IOUT, '(a,2X,E12.4)') "# volumetric_mean(|z|)= ", pdf_all_sum/volume
    write(IOUT, '(a)') "# <z0> <z1> <percentage>"

    do i = 1, nbin
      write(IOUT,'(E12.4, E12.4, E12.4)') z(i-1), z(i), pdf_all(i)
    enddo

    close(IOUT)

  endif

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
