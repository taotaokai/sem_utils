subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_get_dmodel_steepest_descent"
  print '(a)', "    - calculate model udpate direction by steepest descent"
  print '(a)', "      with preconditioning (e.g. mask)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_get_dmodel_steepest_descent \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <model_names> "
  print '(a)', "    <scale_factor> <use_mask> <mask_dir> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  dmodel = scale_factor * [mask] * kernel"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory of proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir:  directory of kernel files "
  print '(a)', "  (string) model_names: comma delimited string, e.g. mu,lamda,rho"
  print '(a)', "  (float) scale_factor:  value to scale kernel values"
  print '(a)', "  (logical) use_mask:  flag if kernel masking is used "
  print '(a)', "  (string) mask_dir:  directory of mask files "
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel "
  print '(a)', "  2. kernel files are <kernel_dir>/proc*_reg1_<model_names>_kernel.bin"
  print '(a)', "  3. mask files are <mask_dir>/proc*_reg1_mask.bin"
  print '(a)', "  4. output files are <out_dir>/proc*_reg1_<model_names>_dmodel.bin"
end subroutine


program xsem_get_dmodel_steepest_descent

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
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: model_names_csv
  real(dp) :: scale_factor
  logical :: use_mask
  character(len=MAX_STRING_LEN) :: mask_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mpi
  integer :: nrank, myrank
  ! model names 
  integer :: nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! mask/dmodel
  real(dp), dimension(:,:,:,:), allocatable :: mask, dmodel

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] xsem_get_dmodel_steepest_descent: check your input arguments."
      call abort_mpi()
    endif
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') kernel_dir 
  read(args(4),'(a)') model_names_csv
  read(args(5),*) scale_factor
  select case (args(6))
    case ('0')
      use_mask = .false.
    case ('1')
      use_mask = .true.
    case default
      if (myrank == 0) then
        print *, '[ERROR]: use_mask must be 0 or 1'
        call abort_mpi() 
      endif
  end select
  read(args(7),'(a)') mask_dir
  read(args(8),'(a)') out_dir

  call synchronize_all()

  !===== parse model_names

  call sem_utils_delimit_string(model_names_csv, ',', model_names, nmodel)

  if (myrank == 0) then
    print *, '# nmodel=', nmodel
    print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel)
  endif

  !===== calculate model update diretion from kernel (steepest descent)

  ! get mesh info: nspec
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  call synchronize_all()

  ! intialize arrays 
  allocate(dmodel(NGLLX,NGLLY,NGLLZ,NSPEC))
  if (use_mask) then
    allocate(mask(NGLLX,NGLLY,NGLLZ,NSPEC))
  endif

  do iproc = myrank, (nproc-1), nrank

    print *, '# iproc=', iproc

    ! read kernels
    do i = 1, nmodel

      ! steepest descent: dmodel is propotional to kernel
      call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, &
        trim(model_names(i))//"_kernel", dmodel)

      ! apply kernel mask
      if (use_mask) then
        call sem_io_read_gll_file_1(mask_dir, iproc, iregion, 'mask', mask)
        dmodel = mask * dmodel
      endif

      ! apply scale factor
      dmodel = scale_factor * dmodel
  
      ! write
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        trim(model_names(i))//"_dmodel", dmodel)

      enddo ! do i

  enddo ! do iproc

  !====== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
