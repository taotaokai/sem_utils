subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_vti_3pars_scale_max_dlnvs"
  print '(a)', "    - scale vti (precondioned-)kernel (vp2, vsv2, vsh2)"
  print '(a)', "      to have some specified maximum Vs_voigt perturbation"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_vti_3pars_scale_max_dlnvs"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <kernel_suffix> "
  print '(a)', "    <max_dlnvs> <out_dir> <out_suffix>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_[vp2,vsv2,vsh2]_kernel.bin"
  print '(a)', "  (string) kernel_suffix:  proc*_reg1_[vp2,vsv2,vsh2]<kernel_suffix>.bin"
  print '(a)', "  (float) max_dlnvs:  maximum shear velocity perturbation to scale to"
  print '(a)', "  (string) out_dir:  output directory for new model"
  print '(a)', "  (string) out_suffix:  output file names proc*_reg1_[vp2,vsv2,vsh2]<out_suffix>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. thomsen's anisotropic parameters: eps = delta = 0 => C11 = C33, C13 = C11 - 2*C44"
  print '(a)', "     vp2_kernel -> C11,C33, vsv2_kernel -> C44, vsh2_kernel -> C66"

end subroutine


program xsem_kernel_vti_3pars_scale_max_dlnvs

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
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: kernel_suffix
  real(dp) :: max_dlnvs_to_scale
  character(len=MAX_STRING_LEN) :: out_dir
  character(len=MAX_STRING_LEN) :: out_suffix

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vsv2, vsh2
  ! kernel
  real(dp), dimension(:,:,:,:), allocatable :: vp2_kernel, vsv2_kernel, vsh2_kernel
  ! scale
  real(dp) :: max_dlnvs, max_dlnvs_all
  real(dp) :: scale_factor

  !===== start MPI
  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      print *, "[ERROR] check your inputs."
      call abort_mpi()
    endif
  endif

  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),*) nproc
  read(args(2),'(a)') mesh_dir 
  read(args(3),'(a)') model_dir 
  read(args(4),'(a)') kernel_dir 
  read(args(5),'(a)') kernel_suffix
  read(args(6),*) max_dlnvs_to_scale
  read(args(7),'(a)') out_dir
  read(args(8),'(a)') out_suffix

  call synchronize_all()

  !===== Get step length first

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(vsv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vp2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  ! get maximum velocity perturbation from dmodel
  max_dlnvs = 0.0_dp

  call synchronize_all()

  do iproc = myrank, (nproc-1), nrank
    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh2)
    ! convert to squared values
    vsv2 = vsv2**2
    vsh2 = vsh2**2
    ! read kernels
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsv2'//trim(kernel_suffix), vsv2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsh2'//trim(kernel_suffix), vsh2_kernel)
    ! calculate maximum absolute relative perturbation of velocities, only use vsv and vsh
    max_dlnvs = max(max_dlnvs, 0.5*maxval(abs(vsv2_kernel)/vsv2), 0.5*maxval(abs(vsh2_kernel)/vsh2))
  enddo ! do iproc

  call synchronize_all()
  call max_all_dp(max_dlnvs, max_dlnvs_all)

  ! determine scale factor 
  if (myrank == 0) then
    scale_factor = max_dlnvs_to_scale / max_dlnvs_all
  endif

  call synchronize_all()

  call bcast_all_singledp(scale_factor)
  print *, "scale_factor = ", scale_factor

  !====== create new model
  call synchronize_all()

  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read kernels
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vp2'//trim(kernel_suffix), vp2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsv2'//trim(kernel_suffix), vsv2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsh2'//trim(kernel_suffix), vsh2_kernel)

    ! write scaled kernels 
    call sem_io_write_gll_file_1(out_dir, iproc, iregion,  'vp2'//trim(out_suffix), scale_factor*vp2_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv2'//trim(out_suffix), scale_factor*vsv2_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh2'//trim(out_suffix), scale_factor*vsh2_kernel)

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
