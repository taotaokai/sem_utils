subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_create_model_homogeneous"
  print '(a)', "    - create a homogeneous SEM model " 
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_create_model_homogeneous \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_name> <model_value> "
  print '(a)', "    <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  The model created is homogeneous with a constant value"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_name: one model name, e.g. vs "
  print '(a)', "  (float) model_value:  the constant model value"
  print '(a)', "  (string) out_dir:  output directory of new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_add_model_horizontal_layer
  
  use sem_constants
  use sem_io
  use sem_mesh

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 6
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_name, out_dir
  integer :: nproc
  real(dp) :: model_value

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), allocatable :: model_gll(:,:,:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_horizontal_layer: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_name
  read(args(5), *) model_value
  read(args(6), '(a)') out_dir

  !===== loop each mesh slice
  do iproc = 0, (nproc - 1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! set model gll
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(NGLLX,NGLLY,NGLLZ,nspec))
    model_gll = model_value

    ! write out model file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, &
      model_name, model_gll)

  enddo ! iproc

end program
