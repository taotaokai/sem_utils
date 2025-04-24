subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_threshold_volume_percentage"
  print '(a)', "    - threshold gll file acccording to volume distribution of amplitudes"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_statis \"
  print '(a)', "    <nproc> <mesh_dir> <hess_dir> <hess_name> "
  print '(a)', "    <nbin> <threshold_percentage> "
  print '(a)', "    <out_dir> <out_name>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) hess_dir:  directory holds proc*_reg1_<hess_name>.bin"
  print '(a)', "  (string) hess_name:  e.g. hess_diag "
  print '(a)', "  (int) nbin:  number of bins"
  print '(a)', "  (float) threshold_percentage: integral(|z|dV, |z|>=z_threshold)/int(|z|dV) <= threshold_percentage, like 0.95"
  print '(a)', "  (string) out_dir:  output directory for threshold inverse hess diag files "
  print '(a)', "  (string) out_name: output gll file as proc*_reg1_<out_name>.bin "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. output is a list of <z0> <z1> <percentage_of_volumetric_amplitude_integral>"
  print '(a)', "  3. hess value |z| <= z_threshold is water leveled to z_threshold"

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_threshold_volume_percentage
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: hess_dir, hess_name
  integer :: nbin
  real(dp) :: threshold_percentage 
  character(len=MAX_STRING_LEN) :: out_dir, out_name

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  real(dp), allocatable :: gll_volume(:,:,:,:)
  ! hess_diag
  real(dp), allocatable :: hess(:,:,:,:)
  ! PDF
  real(dp) :: zmin, zmax, min_i, max_i, bin_size
  real(dp) :: zmin_all, zmax_all
  real(dp), allocatable :: z(:), pdf(:), pdf_all(:)
  real(dp) :: pdf_all_sum, volume, volume_sum
  ! threshold
  real(dp) :: z_threshold

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
  read(args(3), '(a)') hess_dir
  read(args(4), '(a)') hess_name
  read(args(5), *) nbin
  read(args(6), *) threshold_percentage
  read(args(7), '(a)') out_dir
  read(args(8), '(a)') out_name

  ! validate input
  if (myrank == 0) then
    if (nbin <= 1) then
      print *, "[ERROR] nbin must be greater than 1."
      call abort_mpi()
    endif
    if (threshold_percentage<=0.0_dp .or. threshold_percentage>1.0_dp) then
      print *, "[ERROR] threshold_percentage must be between 0 and 1."
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  !===== loop each mesh/hess slice

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! initialize arrays
  allocate(hess(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gll_volume(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(z(0:nbin), pdf(nbin))
  if (myrank == 0) allocate(pdf_all(nbin))

  ! get hess value range
  zmax = -1.0 * huge(1.0_dp)
  zmin = huge(1.0_dp)
  do iproc = myrank, (nproc-1), nrank
    ! read hess
    call sem_io_read_gll_file_1(hess_dir, iproc, iregion, hess_name, hess)
    hess = abs(hess)
    max_i = maxval(hess)
    min_i = minval(hess)
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

  ! maximum amplitude
  !zmax = max(abs(zmax_all), abs(zmin_all))
  ! mesh of bins : (0, max|z|)
  bin_size = (zmax_all - zmin_all) / nbin
  z = (/(i*bin_size, i=0,nbin)/) + zmin_all

  pdf = 0.0_dp
  volume = 0.0_dp
  do iproc = myrank, (nproc-1), nrank
    ! read hess 
    call sem_io_read_gll_file_1(hess_dir, iproc, iregion, hess_name, hess)
    ! get gll_volume
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_mesh_gll_volume(mesh_data, gll_volume)
    ! get volumetric integral of hess amplitudes within each bin 
    hess = abs(hess)
    do i = 1, nbin
      pdf(i) = pdf(i) + sum(gll_volume*hess, mask=(hess>=z(i-1) .and. hess<z(i)) )
    enddo
    ! get volume
    volume = volume + sum(gll_volume)
  enddo
  
  call synchronize_all()

  call sum_all_dp_n(pdf, pdf_all, nbin)
  call sum_all_dp(volume, volume_sum)

  call synchronize_all()

  !------ get threshold amplitude
  if (myrank == 0) then
    ! normalize PDF
    pdf_all_sum = sum(pdf_all)
    pdf_all = pdf_all / pdf_all_sum
    ! get threshold amplitude
    do i = 1,nbin
      if (sum(pdf_all(1:i)) >= (1.0-threshold_percentage)) then
        exit
      endif
    enddo
    z_threshold = z(i)
  endif

  call bcast_all_singledp(z_threshold)

  call synchronize_all()

  !------ write out resutls: pdf
  if (myrank == 0) then
    open(IOUT, file=trim(out_dir)//"/xsem_threshold_percentage.log", status='unknown', &
      form='formatted', action='write', iostat=ier)

    if (ier /= 0) then
      write(*,*) '[ERROR] failed to open log file ', trim(out_dir)//"/xsem_threshold_percentage.log"
      call abort_mpi() 
    endif

    write(IOUT, '(a)') "# PDF of gll hess"
    write(IOUT, '(a,2X,I4)') "# nproc= ", nproc
    write(IOUT, '(a,2X,a)') "# mesh_dir= ", trim(mesh_dir)
    write(IOUT, '(a,2X,a)') "# hess_dir= ", trim(hess_dir)
    write(IOUT, '(a,2X,a)') "# hess_name= ", trim(hess_name)
    write(IOUT, '(a,2X,I4)') "# nbin= ", nbin
    write(IOUT, '(a,2X,E12.4)') "# threshold_percentage = ", threshold_percentage 
    write(IOUT, '(a,2X,E12.4)') "# z_threshold = ", z_threshold
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

  !------ water level hess
  print *, 'z_threshold = ', z_threshold
  do iproc = myrank, (nproc-1), nrank
    print *, 'iproc = ', iproc
    ! read hess 
    call sem_io_read_gll_file_1(hess_dir, iproc, iregion, hess_name, hess)
    ! water level
    hess = abs(hess)
    where (hess<z_threshold) hess = z_threshold
    ! write out inverse water leveled hess diag
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, trim(out_name), 1.0/hess)
  enddo

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
