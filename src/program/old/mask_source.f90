! make mask file for event kernel, e.g. mask source / shallow mantle

subroutine selfdoc()
  print *, 'Usage: xmask_kernel <width(km)> <source_vtk> <topo_dir> <out_dir> <mask_name>'
  print *, 'files read: <kernel_dir>/proc000***_reg1_solver_data.bin'
  print *, 'files written: <out_dir>/proc000***_reg1_<mask_name>.bin'
  stop
end subroutine

program make_sources_mask

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IIN,IOUT, &
    NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL

  use sem_IO

  implicit none

  !---- parameters 
  integer, parameter :: NGLLCUBE = 125
  real(RK), parameter :: EARTH_R_KM = 6371.0

  !---- declare variables
  integer :: ier

  ! command line arguments
  integer :: nargs
  character(len=MAX_STRING_LEN) :: arg
  real(RK) :: mask_width 
  integer :: nslice
  character(len=MAX_STRING_LEN) :: source_file 
  character(len=MAX_STRING_LEN) :: topo_dir 
  character(len=MAX_STRING_LEN) :: outdir
  
  ! source
  character(len=MAX_STRING_LEN) :: tmp 
  integer :: isrc, nsrc
  real(RK), allocatable :: xyz_src(:,:)
  
  ! mesh topo
  character(len=MAX_STRING_LEN) :: topo_file 
  integer :: islice, ispec, iglob, igll, nspec, nglob, slice_id
  real(RK), allocatable :: xyz_glob(:,:)
  integer, allocatable :: ibool(:,:)

  ! source mask
  real(RK) :: dist_sq, weight
  character(len=MAX_STRING_LEN) :: mask_file 
  real(RK), allocatable :: mask_source(:,:)

!---- get input parameters

  ! check argument number
  nargs = command_argument_count() 
  if (nargs < 5) then
    print *, 'Usage: xmask_source <width(km)> <CMTSOLUTION> <NPROC> <topo_dir> <outdir>'
    stop
  endif

  ! mask width
  call get_command_argument(1, arg)
  read(arg,*,iostat=ier) mask_width
  if (ier /= 0) stop 'Error read mask_width'
  if (mask_width < 1E-5) stop 'mask_width too small'
  
  ! source_file
  call get_command_argument(2, arg)
  source_file = trim(arg)

  ! nslice
  call get_command_argument(3, arg)
  read(arg,*,iostat=ier) nslice 
  if (ier /= 0) stop 'Error read nslice' 

  ! topo_dir
  call get_command_argument(4, arg)
  topo_dir = trim(arg) // '/'

  ! outdir
  call get_command_argument(5, arg)
  outdir = trim(arg) // '/'

  ! check
  print *, 'mask_width(km)=', mask_width, &
           ', source_file=', trim(source_file), &
           ', nslice=', nslice, &
           ', topo_dir=', trim(topo_dir), &
           ', outdir=', trim(outdir)

!---- read source file (source.vtk)
  open(unit=IIN, file=source_file)

  ! skip header lines
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  ! nsrc
  read(IIN,*,iostat=ier) tmp, nsrc, tmp
  if (ier /= 0) stop 'Error read nsrc'
  ! xyz_src
  allocate(xyz_src(3,nsrc))
  do isrc = 1, nsrc
    read(IIN,*,iostat=ier) xyz_src(1,isrc), xyz_src(2,isrc), xyz_src(3,isrc)
    if (ier /= 0) stop 'Error read source location' 
  enddo
  close(IIN)

  print *, 'nsrc=', nsrc
  print *, 'xyz_src=', xyz_src 

!---- loop each mesh slice

  ! non-dimensionalize
  mask_width = mask_width / EARTH_R_KM

  do islice = 1, nslice

    slice_id = islice - 1

    !-- read in mesh topo: nspec, nglob, x/y/z_glob(1:nglob), ibool(ngll,nspec)
    write(topo_file,'(a,a,i6.6,a,i1,a)') trim(topo_dir),'/proc',slice_id,'_reg',IREG,'_solver_data.bin'
    open(unit=IIN, file=topo_file, status='old', action='read', iostat=ier, form='unformatted')
    if (ier /= 0) then
      print *,'Error open topo file: ', topo_file
      stop
    endif
    ! read mesh dimension
    read(IIN) nspec
    read(IIN) nglob
    ! read in mesh topology
    if (allocated(xyz_glob)) deallocate(xyz_glob)
    allocate(xyz_glob(3,nglob))
    if (allocated(ibool)) deallocate(ibool)
    allocate(ibool(NGLLCUBE,nspec))
    read(IIN) xyz_glob(1,1:nglob)
    read(IIN) xyz_glob(2,1:nglob)
    read(IIN) xyz_glob(3,1:nglob)
    read(IIN) ibool(1:NGLLCUBE,1:nspec)
    close(IIN)

    ! print
    print *, 'slice_id=', slice_id, 'nspec=',nspec, 'nglob=',nglob 

    !-- loop each element & gll point to set source mask
    if (allocated(mask_source)) deallocate(mask_source)
    allocate(mask_source(NGLLCUBE,nspec))

    do ispec = 1, nspec
      do igll = 1, NGLLCUBE

        iglob = ibool(igll,ispec)

        ! loop each source location to set weight
        weight = 1.0_RK
        do isrc = 1, nsrc
          dist_sq = sum((xyz_glob(:,iglob) - xyz_src(:,isrc))**2)
          weight = weight * (1.0_RK - exp( - (dist_sq/mask_width**2)**2 ))
        enddo

        mask_source(igll,ispec) = weight

      enddo
    enddo

    print *,'min,max(mask_source)=', minval(mask_source), ',', maxval(mask_source)

    !-- save source mask for this slice
    write(mask_file,'(a,a,i6.6,a,i1,a)') trim(outdir),'/proc',slice_id,'_reg',IREG,'_mask_source.bin'
    open(unit=IOUT, file=mask_file, status='unknown', iostat=ier, form='unformatted')
    if (ier /= 0) then
      print *, 'Error open file for write: ', mask_file
      stop
    endif
    write(IOUT) mask_source(:,:)
    close(IOUT)

  enddo ! islice 

!---- clean up

end program make_sources_mask
