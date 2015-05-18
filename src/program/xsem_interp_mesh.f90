subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_interp_mesh - interpolate GLL model from one SEM mesh onto a new mesh"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_interp_mesh <old_mesh_dir> <old_model_dir> <nproc_old> "
  print '(a)', "                   <new_mesh_dir> <new_model_dir> <nproc_new> "
  print '(a)', "                   <region_id> <model_tags> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  interpolate GLL model given on one SEM mesh onto  another SEM mesh"
  print '(a)', "    the GLL interpolation is used, which is the SEM basis function."
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) old_mesh_dir: directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) old_model_dir: directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (int) nproc_old: number of slices of the old mesh"
  print '(a)', "  (string) new_mesh_dir: directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) new_model_dir: directory holds output model files"
  print '(a)', "  (int) nproc_new: number of slices of the new mesh"
  print '(a)', "  (int) region_id: region id of mesh(1: crust_mantle; 2/3: outer/inner ocre)" 
  print '(a)', "  (string) model_tags: comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', ""
  print '(a)', "NOTES"
  print '(a)', ""
  print '(a)', "  This program must run in parallel, e.g. mpirun -n <nproc> ..."

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_interp_mesh
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use sem_parallel

  implicit none

  !===== declare variables

  integer, parameter :: dp = kind(1.0_dp)

  !-- command line args
  integer, parameter :: nargs = 8
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: old_mesh_dir, old_model_dir
  integer :: nproc_old
  character(len=MAX_STRING_LEN) :: new_mesh_dir, new_model_dir
  integer :: nproc_new
  integer :: region_id
  character(len=MAX_STRING_LEN) :: model_tags

  !-- local variables
  integer :: i, iproc, ier, nrank, myrank

  !-- model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- mesh
  type(sem_mesh) :: old_mesh, new_mesh
  integer :: igllx, iglly, igllz, ispec
  real(dp), allocatable :: xyz(:,:)

  !-- sem interpolation
  real(CUSTOM_REAL), allocatable :: model_gll(:,:,:,:,:), model_interp(:,:)
  real(CUSTOM_REAL), allocatable :: uvw(:,:), hlagrange(:,:,:,:)
  real(CUSTOM_REAL), allocatable :: misloc(:), misloc_final(:)
  integer, allocatable :: eid(:), eid_final(:)
  integer, allocatable :: locstat(:), locstat_final(:)

  !===== start MPI

  call init_mpi()
  call world_size(nrank)
  call world_rank(myrank)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    if (myrank == 0) then
      call selfdoc()
      stop "[ERROR] xsem_interp_mesh: check your input arguments."
    endif
  endif
  call synchronize_all()

  do i = 1, nargs
    call get_command_argument(i, args(i), status=ier)
  enddo
  read(args(1), '(a)') old_mesh_dir
  read(args(2), '(a)') old_model_dir
  read(args(3), *) nproc_old
  read(args(4), '(a)') new_model_dir
  read(args(5), '(a)') new_mesh_dir
  read(args(6), *) nproc_new
  read(args(7), *) region_id
  read(args(8), '(a)') model_tags

  !-- read model tags
  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  if (myrank == 0) then
    print '(a)', '# nmodel=', nmodel
    print '(a)', '# model_names=', (trim(model_names(i))//" ", i=1,nmodel)
  endif

  !===== loop each slices of the new mesh

  do iproc_new = myrank, (nproc_new - 1), nrank

    !-- read new mesh slice
    call sem_mesh_read(mesh_new, new_mesh_dir, iproc_new, region_id)
    




  enddo ! iproc_new

  !===== locate xyz in the mesh

  ! initialize arrays for each mesh chunk
  allocate(model_gll(NGLLX,NGLLY,NGLLZ,NSPEC,nmodel), &
           model_interp(nmodel,npoint), &
           uvw(3,npoint), &
           hlagrange(NGLLX,NGLLY,NGLLZ,npoint), &
           misloc(npoint), &
           misloc_final(npoint), &
           elem_ind(npoint), &
           elem_ind_final(npoint), &
           statloc(npoint), &
           statloc_final(npoint))

  ! initialize some variables
  statloc_final = -1
  misloc_final = HUGE(1.0_CUSTOM_REAL)
  elem_ind_final = -1
  model_interp = FILLVALUE

  ! typical element size at surface
  typical_size = real( max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, &
                           ANGULAR_WIDTH_ETA_IN_DEGREES_VAL/ NEX_ETA_VAL) &
                       * DEGREES_TO_RADIANS * R_UNIT_SPHERE, kind=CUSTOM_REAL)

  ! loop each mesh chunk
  call sem_mesh_init(mesh_data)
  do iproc = 0, NPROCTOT_VAL-1

    print *, "# iproc=", iproc

    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    call sem_io_read_gll_modeln(model_dir, iproc, iregion, nmodel, model_names &
                                , model_gll)
    call sem_mesh_locate_xyz(mesh_data, npoint, xyz, uvw, hlagrange, misloc &
                             , elem_ind, statloc)

    print *, "# number of located points: ", count(statloc>=0)

    ! interpolate model only on points located inside an element
    do ipoint = 1, npoint

      ! safety check
      if (statloc_final(ipoint) == 1 .and. statloc(ipoint) == 1 ) then
        print *, "WARN: ipoint=", ipoint
        print *, "====: this point is located inside more than one element!"
        print *, "====: some problem may occur."
        print *, "====: only use the first located element."
        cycle
      endif

      ! for point located inside one element but not happened before 
      ! or closer to one element than located before
      if ( (statloc(ipoint) == 1 .and. statloc_final(ipoint) /= 1 ) .or. &
           (statloc(ipoint) == 0 .and. misloc(ipoint) < misloc_final(ipoint)) ) then

        ! interpolate model
        do imodel = 1, nmodel
          model_interp(imodel,ipoint) = sum( &
            hlagrange(:,:,:,ipoint) * model_gll(:,:,:,elem_ind(ipoint),imodel) )
        enddo

        statloc_final(ipoint) = statloc(ipoint)
        misloc_final(ipoint) = misloc(ipoint)
        elem_ind_final(ipoint) = elem_ind(ipoint)

      endif

    enddo ! ipoint

  enddo ! iproc

  ! convert misloc to relative to typical element size
  where (statloc_final /= -1) misloc_final = misloc_final / typical_size

  !===== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
