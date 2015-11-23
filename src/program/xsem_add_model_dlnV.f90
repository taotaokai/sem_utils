subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_model_dlnV "
  print '(a)', "    - add model relative perturbations into SEM GLL files"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_model_dlnV \"
  print '(a)', "    <mesh_dir> <nproc> "
  print '(a)', "    <model_dir> <model_tags> "
  print '(a)', "    <dlnV_dir> <dlnV_tags> <scale_factor>"
  print '(a)', "    <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  new_model = model * (1 + scale_factor * dlnV)"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vs,vp "
  print '(a)', "  (string) dlnV_dir:  directory holds relative model perturbation (proc*_reg1_<dlnV_tag>.bin)"
  print '(a)', "  (string) dlnV_tags:  comma delimited string, e.g. dvs,dvp "
  print '(a)', "  (float) scale_factor:  scaling factor "
  print '(a)', "  (string) out_dir:  output directory for the new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_add_model_dlnV
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir
  integer :: nproc
  character(len=MAX_STRING_LEN) :: model_dir, model_tags
  character(len=MAX_STRING_LEN) :: dlnV_dir, dlnV_tags
  real(dp) :: scale_factor
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! model names
  integer :: nmodel, ndlnV
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)
  character(len=MAX_STRING_LEN), allocatable :: dlnV_names(:)

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), allocatable :: gll_model(:,:,:,:,:), gll_dlnV(:,:,:,:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_dlnV: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), '(a)') dlnV_dir
  read(args(6), '(a)') dlnV_tags
  read(args(7), *) scale_factor
  read(args(8), '(a)') out_dir 

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)
  call sem_utils_delimit_string(dlnV_tags, ',', dlnV_names, ndlnV)

  if (nmodel /= ndlnV) then
    print *, '[ERROR] nmodel /= ndlnV_2', nmodel, ndlnV
    stop
  endif

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel)

  print *, '# ndlnV=', ndlnV
  print *, '# dlnV_names=', (trim(dlnV_names(i))//"  ", i=1,ndlnV)

  !===== loop each mesh/model slice

  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! read gll model
    if (allocated(gll_model)) then
      deallocate(gll_model, gll_dlnV)
    endif
    allocate(gll_model(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
    allocate(gll_dlnV(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
                                model_names, nmodel, gll_model)

    call sem_io_read_gll_file_n(dlnV_dir, iproc, iregion, &
                                dlnV_names, nmodel, gll_dlnV)

    ! calculate perturbed model
    gll_model = gll_model * (1.0 + scale_factor*gll_dlnV)

    ! write out the perturbed model
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, &
      model_names, nmodel, gll_model)

  enddo ! iproc

end program
