subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_interp_mesh - interpolate GLL model from one SEM mesh onto a new mesh"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_interp_mesh <old_mesh_dir> <old_model_dir> <nproc_old> "
  print '(a)', "                   <new_mesh_dir> <new_model_dir> <nproc_new> "
  print '(a)', "                   <iregion> <model_tags> "
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
  print '(a)', "  (int) iregion: region id of mesh(1: crust_mantle; 2/3: outer/inner ocre)" 
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
  integer :: iregion
  character(len=MAX_STRING_LEN) :: model_tags

  !-- local variables
  integer :: i, iproc, ier, nrank, myrank

  !-- model names
  integer :: imodel, nmodel
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  !-- mesh
  type(sem_mesh) :: mesh_old, mesh_new
  integer :: igllx, iglly, igllz, ispec
  real(dp), allocatable :: xyz(:,:)

  !-- sem interpolation
  real(CUSTOM_REAL), allocatable :: model_gll_old(:,:,:,:,:)
  real(CUSTOM_REAL), allocatable :: model_gll_new(:,:,:,:,:), model_interp(:,:)
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
  read(args(7), *) iregion
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
    call sem_mesh_read(mesh_new, new_mesh_dir, iproc_new, iregion)

    !-- initialize arrays of xyz points
    nspec_new = mesh_new%nspec
    nglob_new = mesh_new%nglob
    ngll_new = NGLLX * NGLLY * NGLLZ * nspec_new

    if (allocated(xyz_new)) then
      deallocate(xyz_new, idoubling_new)
    endif
    allocate(xyz_new(3, ngll_new), idoubling_new(ngll_new))

    do ispec = 1, nspec_new
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            igll = igllx + &
                   NGLLX * ( (iglly-1) + &
                   NGLLY * ( (igllz-1) + &
                   NGLLZ * ( (ispec-1))))
            iglob = mesh_new%ibool(igllx, iglly, igllz, ispec)
            xyz_new(:, igll) = mesh_new%xyz(:, iglob)
            idoubling_new(igll) = mesh_new%idoubling(ispec)
          enddo
        enddo
      enddo
    enddo

    !-- initialize arrays
    if (allocated(uvw)) then
      deallocate(uvw, hlagrange, misloc, eid, locstat, model_interp, &
        model_gll_new)
    endif
    allocate(uvw(3,ngll_new), &
             hlagrange(NGLLX,NGLLY,NGLLZ,ngll_new), &
             misloc(ngll_new), &
             eid(ngll_new), &
             locstat(ngll_new), &
             model_interp(nmodel,ngll_new), &
             model_gll_new(nmodel,NGLLX,NGLLY,NGLLZ,nspec_new))

    if (allocated(uvw_final)) then
      deallocate(uvw_final, hlagrange_final, misloc_final, &
                 eid_final, locstat_final)
    endif
    allocate(uvw_final(3,ngll_new), &
             hlagrange_final(NGLLX,NGLLY,NGLLZ,ngll_new), &
             misloc_final(ngll_new), &
             eid_final(ngll_new), &
             locstat_final(ngll_new))

    hlagrange_final = HUGE(1.0_dp)
    misloc_final = HUGE(1.0_dp)
    eid_final = -1
    locstat_final = -1
    
    !-- loop each slices of the old mesh
    do iproc_old = 0, (nproc_old - 1)

      ! read old mesh slice
      call sem_mesh_read(mesh_old, old_mesh_dir, iproc_old, iregion)

      dims_old = [NGLLX, NGLLY, NGLLZ, mesh_old%nspec]

      ! read old model
      if (allocated(model_gll_old)) then
        deallocate(model_gll_old)
      endif
      allocate(model_gll_old(nmodel, NGLLX, NGLLY, NGLLZ, mesh_old%nspec))

      call sem_io_read_gll_n(old_model_dir, iproc, iregion, &
                             nmodel, model_names, dims_old, model_gll_old)

      call sem_mesh_locate_xyz(mesh_old, ngll_new, xyz_new, idoubling_new, &
                               uvw, hlagrange, misloc, eid, locstat)

      ! interpolate model only on points located inside an element
      do igll = 1, ngll_new

        ! safety check
        if (locstat_final(igll) == 1 .and. locstat(igll) == 1 ) then
          print *, "[WARN] igll=", igll
          print *, "------ this point is located inside more than one element!"
          print *, "------ some problem may occur."
          print *, "------ only use the first located element."
          cycle
        endif

        ! for point located inside one element for the first time 
        ! or closer to one element than located before
        if ( (locstat(igll) == 1 .and. locstat_final(igll) /= 1 ) .or. &
             (locstat(igll) == 0 .and. misloc(igll) < misloc_final(igll)) ) &
        then

          ! interpolate model
          do imodel = 1, nmodel
            model_interp(imodel,igll) = sum(hlagrange(:,:,:,igll) * &
              model_gll_old(imodel,:,:,:,eid(igll)))
          enddo

          locstat_final(igll) = locstat(igll)
          misloc_final(igll) = misloc(igll)
          eid_final(igll) = eid(igll)

        endif

      enddo ! igll = 1, ngll_new

    enddo ! iproc_old

    !-- output location information

    ! convert misloc to relative to typical element size
    where (statloc_final /= -1) 
      misloc_final = misloc_final / typical_size

    !-- write out gll files for new mesh slice

    do ispec = 1, nspec_new
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
            igll = igllx + &
                   NGLLX * ( (iglly-1) + & 
                   NGLLY * ( (igllz-1) + & 
                   NGLLZ * ( (ispec-1))))
            if (locstat_final /= -1) then
              model_gll_new(:, igllx, iglly, igllz, ispec) = &
                model_interp(:, igll)
            endif
          enddo
        enddo
      enddo
    enddo

    call sem_io_write_gll_file_n(new_model_dir, iproc_new, iregion, nmodel, &
      model_names, model_gll_new)

  enddo ! iproc_new

  !===== exit MPI
  call synchronize_all()
  call finalize_mpi()

end program
