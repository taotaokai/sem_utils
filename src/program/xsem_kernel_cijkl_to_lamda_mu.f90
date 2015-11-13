subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_cijkl_kernel_to_lamda_mu "
  print '(a)', "    - reduce cijkl kernel to (lamda,mu) kernel"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_cijkl_kernel_to_lamda_mu \"
  print '(a)', "    <nproc> <mesh_dir> <kernel_dir> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  "
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) kernel_dir:  directory holds proc*_reg1_cijkl_kernel.bin"
  print '(a)', "  (string) out_dir:  output directory for lamda,mu_kerenl"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_cijkl_kernel_to_lamda_mu

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables
  ! command line args
  integer, parameter :: nargs = 4
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: kernel_dir
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec

  ! kernel gll 
  real(dp), allocatable :: cijkl_kernel(:,:,:,:,:)
  real(dp), allocatable :: lamda_kernel(:,:,:,:), mu_kernel(:,:,:,:)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_cijkl_kernel_to_lamda_mu: check your inputs."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), *) mesh_dir
  read(args(3), *) kernel_dir
  read(args(4), *) out_dir 

  !====== loop model slices 
  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    nspec = mesh_data%nspec

    ! read gll model
    if (allocated(cijkl_kernel)) then
      deallocate(cijkl_kernel, lamda_kernel, mu_kernel)
    endif
    allocate(cijkl_kernel(21,NGLLX,NGLLY,NGLLZ,nspec))
    allocate(lamda_kernel(NGLLX,NGLLY,NGLLZ,nspec))
    allocate(mu_kernel(NGLLX,NGLLY,NGLLZ,nspec))

    ! read kernel gll
    call sem_io_read_gll_file_cijkl(kernel_dir, iproc, iregion, cijkl_kernel)

    ! reduce cijkl kernel to lamda,mu kernel 
    lamda_kernel = cijkl_kernel(1,:,:,:,:) + &
                   cijkl_kernel(2,:,:,:,:) + &
                   cijkl_kernel(3,:,:,:,:) + &
                   cijkl_kernel(7,:,:,:,:) + &
                   cijkl_kernel(8,:,:,:,:) + &
                   cijkl_kernel(12,:,:,:,:)

    mu_kernel = cijkl_kernel(1,:,:,:,:)  + &
                cijkl_kernel(7,:,:,:,:)  + &
                cijkl_kernel(12,:,:,:,:) + &
                cijkl_kernel(16,:,:,:,:) + &
                cijkl_kernel(19,:,:,:,:) + &
                cijkl_kernel(21,:,:,:,:)

    ! write out lamda,mu kernel
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        'lamda_kernel', lamda_kernel)
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
        'mu_kernel', mu_kernel)

  enddo ! iproc

end program
