subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_inner_product "
  print '(a)', "    - compute inner product between two SEM gll models"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_inner_product \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir1> <model_dir2> <model_tags> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  the inner product are defined as the volume integral of "
  print '(a)', "     the product between the gll model 1 and 2 "
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir1/2:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. the unit of volume is km^3 "
  print '(a)', "  2. do NOT care about the unit of model files, "
  print '(a)', "    , only values are used. "
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_vertical_slice
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir1, model_dir2, model_tags
  integer :: nproc

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec
  ! model
  real(dp), allocatable :: gll_volume(:,:,:,:)
  real(dp), allocatable :: gll_model_1(:,:,:,:,:), gll_model_2(:,:,:,:,:)

  !-- inner product
  real(dp), allocatable :: inner_product(:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_plume: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir1
  read(args(4), '(a)') model_dir2
  read(args(5), '(a)') model_tags

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel) 

  !===== loop each mesh/model slice

  ! initialize arrays
  allocate(inner_product(nmodel,0:(nproc-1)))

  do iproc = 0, (nproc - 1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! calculate gll volume
    if (allocated(gll_volume)) then
      deallocate(gll_volume)
    endif
    allocate(gll_volume(NGLLX,NGLLY,NGLLZ,nspec))
   
    call sem_mesh_gll_volume(mesh_data, gll_volume)

    ! read gll model
    if (allocated(gll_model_1)) then
      deallocate(gll_model_1, gll_model_2)
    endif
    allocate(gll_model_1(nmodel,NGLLX,NGLLY,NGLLZ,nspec))
    allocate(gll_model_2(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

    call sem_io_read_gll_file_n(model_dir1, iproc, iregion, &
                                model_names, nmodel, gll_model_1)

    call sem_io_read_gll_file_n(model_dir2, iproc, iregion, &
                                model_names, nmodel, gll_model_2)

    ! calculate inner product for each model 
    do imodel = 1, nmodel
      inner_product(imodel,iproc) = sum( gll_volume * &
        gll_model_1(imodel,:,:,:,:) * gll_model_1(imodel,:,:,:,:))
    enddo

  enddo ! iproc

  !===== print out results

  do imodel = 1, nmodel
    print "(a,' = ',E15.7)", trim(model_names(imodel)), &
      sum(inner_product(imodel,:)) * R_EARTH_KM**3
  enddo

  print "('sum = ',E15.7)", sum(inner_product) * R_EARTH_KM**3

end program
