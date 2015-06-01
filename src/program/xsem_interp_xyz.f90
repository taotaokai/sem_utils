subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_interp_xyz - interpolate SEM model on given points "
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_interp_xyz \ "
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <xyz_list> <out_list>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  interpolate SEM mesh and associated model values onto the given xyz list"
  print '(a)', "    the SEM GLL interpolation function is used."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices "
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (string) xyz_list:  input file name, "
  print '(a)', "    each row: <x> <y> <z>"
  print '(a)', "  (string) out_list:  output file name, "
  print '(a)', "    each row: <x> <y> <z> misloc <misloc> [<tag> <val>]... "
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_interp_xyz
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !---- declare variables
  ! command line args
  integer, parameter :: nargs = 5
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, xyz_list, out_list

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier 

  integer :: imodel, nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  integer :: ipoint, npoint
  real(CUSTOM_REAL), allocatable :: xyz(:,:)
  character(len=MAX_STRING_LEN), allocatable :: xyz_strings(:) 

  type (mesh) :: mesh_data
  real(CUSTOM_REAL), allocatable :: model_gll(:,:,:,:,:), model_interp(:,:)
  real(CUSTOM_REAL), allocatable :: uvw(:,:), hlagrange(:,:,:,:), misloc(:)
  integer, allocatable :: elem_ind(:), stat_loc(:), stat_loc_final(:)

  !---- read command line arguments
  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
    if (i <= nargs .and. len_trim(args(i)) == 0) then
      call selfdoc()
      stop
    endif
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), '(a)') model_dir
  read(args(3), '(a)') model_tags
  read(args(4), '(a)') xyz_list
  read(args(5), '(a)') out_list

  !---- read xyz list
  call sem_utils_read_line(xyz_list, xyz_strings, npoint)

  allocate(xyz(3, npoint))

  do i = 1, npoint
    read(xyz_strings(i),*,iostat=ier) xyz(1,i), xyz(2,i), xyz(3,i)
    if (ier /= 0) then
      print '(a)', "ERROR: read list at line ", i
      stop
    endif
  enddo

  print '(a)', '# npoint=', npoint
 
  !---- read model tags 
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print '(a)', '# nmodel=', nmodel
  print '(a)', '# model_names=', model_names

  !---- locate xyz in the mesh 
  call sem_constants_set(iregion)

  ! initialize arrays for each mesh chunk
  allocate(model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel), &
           model_interp(nmodel,npoint), &
           uvw(NDIM,npoint), &
           hlagrange(NGLLX,NGLLY,NGLLZ,npoint), &
           misloc(npoint), &
           elem_ind(npoint), &
           stat_loc(npoint), &
           stat_loc_final(npoint))

  stat_loc_final = -1
  model_interp = 0.0_CUSTOM_REAL

  ! loop each mesh chunk
  do iproc = 0, NPROCTOT_VAL-1
  !do iproc = 0, 0

    print '(a)', '# iproc=', iproc

    call sem_mesh_init(mesh_data)
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_io_read_gll_modeln(model_dir, iproc, iregion, nmodel, model_names &
                                , model_gll)
    call sem_mesh_locate_xyz(mesh_data, npoint, xyz, uvw, hlagrange, misloc &
                             , elem_ind, stat_loc)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! for point located inside an element
      if (stat_loc(ipoint) == 1) then

        ! should not happen unless exactly on the element surface
        if (stat_loc_final(ipoint) == 1) then
          print '(a)', "WARN: more than one elements locate point ", ipoint
          print '(a)', "WARN: only use the first located element."
          cycle
        endif

        ! interpolate model
        do imodel = 1, nmodel
          model_interp(imodel,ipoint) = sum( &
            hlagrange(:,:,:,ipoint) * model_gll(:,:,:,elem_ind(ipoint),imodel) )
        enddo

        stat_loc_final(ipoint) = 1

      endif

    enddo ! ipoint

  enddo ! iproc

  !---- write out results
  open(unit=IOUT, file=out_list, status='unknown', iostat=ier)
  if (ier /= 0) then
    print '(a)', 'ERROR: failed to open file: ', out_list
    stop
  endif

  do ipoint = 1, npoint
    write(IOUT, '(3E15.7,a,E15.7)', advance='no') &
      xyz(:,ipoint),' misloc ', misloc(ipoint)
    do imodel = 1, nmodel
      write(IOUT, '(2x,a,2x,E15.7)', advance='no') &
        trim(model_names(i)), model_interp(imodel, ipoint)
    enddo
    write(IOUT,'(a)') ""
  enddo
  close(IOUT)

end program
