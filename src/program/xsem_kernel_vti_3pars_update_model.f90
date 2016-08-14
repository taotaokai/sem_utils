subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_kernel_vti_3pars_update_model"
  print '(a)', "    - add (scaled) vti (precondioned-)kernel (vp2, vsv2, vsh2)"
  print '(a)', "      to GLL vti model file (vpv,vph,vsv,vsh,eta,rho)"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_kernel_vti_3pars_update_model"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <kernel_dir> <kernel_suffix>"
  print '(a)', "    <scale_vp2> <scale_vsv2> <scale_vsh2> <dlnv_limit> <output_dmodel> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  vph2,vpv2 += scale_vp2 * vp2_kernel"
  print '(a)', "  vsv2,vsh2 += scale_vsv2,vsh2 * vsv2,vsh2_kernel"
  print '(a)', "  vf2 +=  scale_vp2*vp2_kernel - 2*scale_vsv2*vsv2_kernel"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_[vp2,vsv2,vsh2]_kernel.bin"
  print '(a)', "  (string) kernel_suffix:  proc*_reg1_[vp2,vsv2,vsh2]<kernel_suffix>.bin"
  print '(a)', "  (float) scale_vp2:  scale factor for vp2_kernel"
  print '(a)', "  (float) scale_vsv2:  scale factor for vsv2_kernel"
  print '(a)', "  (float) scale_vsh2:  scale factor for vsh2_kernel"
  print '(a)', "  (float) dlnv_limit:  limit of maximum velocity change (no limit if < 0)"
  print '(a)', "  (int) output_dmodel(0/1): =1 write out [vp2,vsv2,vsh2]_dmodel.bin, =0 do not"
  print '(a)', "  (string) out_dir:  output directory for new model"
  print '(a)', ""
  print '(a)', "note"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"
  print '(a)', "  2. thomsen's anisotropic parameters: eps = delta = 0 => C11 = C33, C13 = C11 - 2*C44"
  print '(a)', "     vp2_kernel -> C11,C33, vsv2_kernel -> C44, vsh2_kernel -> C66"
  print '(a)', "  3. for isotropic element the program will reduce the vti kernel to isotropic kernel and alos enforce the output be isotropic."

END SUBROUTINE


program xsem_kernel_vti_3pars_update_model

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 11
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: kernel_suffix
  real(dp) :: scale_vp2
  real(dp) :: scale_vsv2
  real(dp) :: scale_vsh2
  real(dp) :: dlnv_limit
  integer :: output_dmodel 
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mpi
  integer :: myrank, nrank
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: ispec, nspec
  ! model
  real(dp), dimension(:,:,:,:), allocatable :: vph2, vpv2, vsv2, vsh2, vf2, rho
  ! kernel
  real(dp), dimension(:,:,:,:), allocatable :: vp2_kernel, vsv2_kernel, vsh2_kernel
  ! scale rho 
  real(dp), parameter :: scale_rho = 0.33
  real(dp), dimension(:,:,:,:), allocatable :: dln_rho
  ! dlnv_limit
  real(dp), dimension(:,:,:,:), allocatable :: dlnv, ones

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
  read(args(6),*) scale_vp2
  read(args(7),*) scale_vsv2
  read(args(8),*) scale_vsh2
  read(args(9),*) dlnv_limit
  read(args(10),*) output_dmodel
  read(args(11),'(a)') out_dir

  call synchronize_all()

  ! get mesh geometry
  if (myrank == 0) then
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)
    nspec = mesh_data%nspec
  endif
  call bcast_all_singlei(nspec)

  ! intialize arrays 
  allocate(vph2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vf2(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(vp2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsv2_kernel(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vsh2_kernel(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dln_rho(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(dlnv(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ones(NGLLX,NGLLY,NGLLZ,nspec))
 
  !====== create new model
  do iproc = myrank, (nproc-1), nrank
  
    print *, "iproc = ", iproc

    ! read in mesh 
    call sem_mesh_read(mesh_dir, myrank, iregion, mesh_data)

    ! read models
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vph', vph2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vpv', vpv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsv', vsv2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'vsh', vsh2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'eta', vf2)
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, 'rho', rho)
    ! convert to squared value and also vf2
    vph2 = vph2**2              ! = C11/rho
    vpv2 = vpv2**2              ! = C33/rho
    vsv2 = vsv2**2              ! = C44/rho
    vsh2 = vsh2**2              ! = C66/rho
    vf2 = vf2*(vph2 - 2.0*vsv2) ! = C13/rho

    ! read kernels
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion,  'vp2'//kernel_suffix, vp2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsv2'//kernel_suffix, vsv2_kernel)
    call sem_io_read_gll_file_1(kernel_dir, iproc, iregion, 'vsh2'//kernel_suffix, vsh2_kernel)

    ! scale kernels
    vp2_kernel = scale_vp2 * vp2_kernel
    vsv2_kernel = scale_vsv2 * vsv2_kernel
    vsh2_kernel = scale_vsh2 * vsh2_kernel

    ! limit maximum dlnv
    if (dlnv_limit > 0.0) then
      ones = 1.0_dp

      dlnv = 0.5*abs(vp2_kernel/vph2)
      where (dlnv > dlnv_limit) vp2_kernel = 2.0*dlnv_limit*vph2*sign(ones,vp2_kernel)

      dlnv = 0.5*abs(vp2_kernel/vpv2)
      where (dlnv > dlnv_limit) vp2_kernel = 2.0*dlnv_limit*vpv2*sign(ones,vp2_kernel)

      dlnv = 0.5*abs(vsv2_kernel/vsv2)
      where (dlnv > dlnv_limit) vsv2_kernel = 2.0*dlnv_limit*vsv2*sign(ones,vsv2_kernel)

      dlnv = 0.5*abs(vsh2_kernel/vsh2)
      where (dlnv > dlnv_limit) vsh2_kernel = 2.0*dlnv_limit*vsh2*sign(ones,vsh2_kernel)
    endif

    ! check if element is isotropic
    do ispec = 1, nspec
      if (.not. mesh_data%ispec_is_tiso(ispec)) then
        ! enforce model to be isotropic
        vph2(:,:,:,ispec) = vpv2(:,:,:,ispec)
        vsh2(:,:,:,ispec) = vsv2(:,:,:,ispec)
        vf2(:,:,:,ispec) = vpv2(:,:,:,ispec) - 2.0*vsv2(:,:,:,ispec)
        ! reduce vti kernel to iso kernel
        vsv2_kernel(:,:,:,ispec) = vsv2_kernel(:,:,:,ispec) + vsh2_kernel(:,:,:,ispec)
        vsh2_kernel(:,:,:,ispec) = vsv2_kernel(:,:,:,ispec)
      endif
    enddo

    ! scale
    ! voigt_G = (2*C44 + C66)/3
    dln_rho = scale_rho*(2.0*vsv2_kernel + vsh2_kernel)/(2.0*vsv2 + vsh2)/2.0
    rho = rho*(1.0 + dln_rho)
    print *, "dln_rho: min/max = ", minval(dln_rho), maxval(dln_rho)

    vph2 = vph2 + vp2_kernel
    vpv2 = vpv2 + vp2_kernel
    vsv2 = vsv2 + vsv2_kernel
    vsh2 = vsh2 + vsh2_kernel
    vf2  = vf2  + (vp2_kernel - 2.0*vsv2_kernel)

    ! write new models (for model_gll)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vph', sqrt(vph2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vpv', sqrt(vpv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv', sqrt(vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh', sqrt(vsh2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'eta', vf2/(vph2 - 2.0*vsv2))
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'rho', rho)

    ! write out dmodel
    if (output_dmodel == 1) then
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vp2_dmodel', vp2_kernel)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsv2_dmodel', vsv2_kernel)
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, 'vsh2_dmodel', vsh2_kernel)
    endif

  enddo

  !====== Finalize
  call synchronize_all()
  call finalize_mpi()

end program
