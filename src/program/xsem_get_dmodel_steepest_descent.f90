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
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <out_dir>"
  print '(a)', "    <model_tags> <scale_factor> <use_mask> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  dmodel = scale_factor * [mask] * kernel"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir:  directory holds proc***_reg1_***_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory"
  print '(a)', "  (string) model_tags: comma delimited string, e.g. mu,lamda,rho"
  print '(a)', "  (float) scale_factor:  factor to scale kernel values"
  print '(a)', "  (logical) use_mask:  flag if kernel masking is used (mask files should locate in kernel_dir)"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. _kernel will be appended to model_tags to get kernel files"
  print '(a)', "  2. _dmodel will be appended to model_tags for output files"
end subroutine


program xsem_get_dmodel_steepest_descent

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 7
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: model_tags
  real(dp) :: scale_factor
  logical :: use_mask

  ! model names
  integer :: nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! mask/dmodel
  real(dp), dimension(:,:,:,:), allocatable :: mask, dmodel

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_get_dmodel_steepest_descent: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') kernel_dir 
  read(args(4),'(a)') out_dir
  read(args(5),'(a)') model_tags
  read(args(6),*) scale_factor 
  select case (args(7))
    case ('0')
      use_mask = .false.
    case ('1')
      use_mask = .true.
    case default
      print *, '[ERROR]: use_mask must be 0 or 1'
      stop
  end select

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel)

  !===== calculate model update diretion from kernel (steepest descent)

  ! get mesh info: nspec
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  ! intialize arrays 
  allocate(dmodel(NGLLX,NGLLY,NGLLZ,NSPEC))
  if (use_mask) then
    allocate(mask(NGLLX,NGLLY,NGLLZ,NSPEC))
  endif

  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! read kernels
    do i = 1, nmodel

      call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, &
        trim(model_names(i))//"_kernel", dmodel)

      ! apply kernel mask
      if (use_mask) then
        call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'mask', mask)
        dmodel = mask * dmodel
      endif

      ! apply scale factor
      dmodel = scale_factor * dmodel
  
      ! write
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        trim(model_names(i))//"_dmodel", dmodel)

      enddo ! do i

  enddo ! do iproc

end program
