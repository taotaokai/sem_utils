subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_lamda_mu_to_dmodel "
  print '(a)', "    - make isotropic model perturbation directions(dlamda, dmu, drho)"
  print '(a)', "      from the isotropic kernels"  
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_lamda_mu_to_dmodel \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <out_dir>"
  print '(a)', "    <model_tags> <scale_factor> <use_kernel_mask> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  dmodel = scale_factor * [kernel_mask] * kernel"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_[lamda,mu,rho]_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for dmodel"
  print '(a)', "  (string) model_tags: comma delimited model names to be used, e.g. mu_kernel,lamda_kernel"
  print '(a)', "  (float) scale_factor:  factor to scale kernel values"
  print '(a)', "  (logical) use_kernel_mask:  flag if kernel mask is used (mask files should locate in kernel_dir)"

  print '(a)', ""
end subroutine


program xsem_kernel_lamda_mu_to_dmodel

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
  logical :: use_kernel_mask

  ! model names
  integer :: nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! mask/kernel/dmodel
  real(dp), dimension(:,:,:,:), allocatable :: mask, dmu, dlamda, drho

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_kernel_lamda_mu_to_dmodel_sd: check your input arguments."
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
      use_kernel_mask = .false.
    case ('1')
      use_kernel_mask = .true.
    case default
      print *, '[ERROR]: use_kernel_mask must be 0 or 1'
      stop
  end select

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel)

  !===== create dmodel from kernel (steepest descent)

  ! get mesh info: nspec
  call sem_mesh_read(mesh_dir, 0, iregion, mesh_data)
  nspec = mesh_data%nspec

  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! intialize mode/dmodel arrays 
    if (.not. allocated(dlamda)) then
      allocate(dmu(NGLLX,NGLLY,NGLLZ,NSPEC), &
            dlamda(NGLLX,NGLLY,NGLLZ,NSPEC), &
              drho(NGLLX,NGLLY,NGLLZ,NSPEC))
      if (use_kernel_mask) then
        allocate(mask(NGLLX,NGLLY,NGLLZ,NSPEC))
      endif
    endif

    ! read kernels
    dlamda = 0.0_dp
    dmu = 0.0_dp
    drho = 0.0_dp
    do i = 1, nmodel
      select case (model_names(i))
        case ('lamda_kernel')
          call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'lamda_kernel', dlamda)
        case ('mu_kernel')
          call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'mu_kernel', dmu)
        case ('rho_kernel')
          call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'rho_kernel', drho)
        case DEFAULT
          print *, '[WARNING] invalid model name: ', model_names(i)
      end select
    enddo

    ! apply kernel mask
    if (use_kernel_mask) then
      call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'mask', mask)
      dlamda = mask * dlamda
      dmu = mask * dmu
      drho = mask * drho
    endif

    ! apply scale factor
    dlamda = scale_factor * dlamda
    dmu = scale_factor * dmu
    drho = scale_factor * drho
  
    ! write dmodel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'lamda_dmodel', dlamda)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'mu_dmodel', dmu)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho_dmodel', drho)

  enddo ! do iproc

end program
