subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_mesh_get_slice_center"
  print '(a)', "    - make GLL files of the slice number"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_mesh_get_slice_center \"
  print '(a)', "    <mesh_dir> <iproc>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
end subroutine


program xsem_mesh_get_slice_center

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 2
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir
  integer :: nproc

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, ier
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: num, id, iproc, ispec, nspec, iglob
  ! center xyz
  real(dp), allocatable :: xyz(:,:)
  real(dp) :: xyz_avg(3), lat, lon

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    stop "[ERROR] check your inputs."
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1),'(a)') mesh_dir 
  read(args(2),*) nproc

  !====== get centers of each slice 
  do iproc = 0, (nproc - 1)

    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    nspec = mesh_data%nspec

    ! separate the 3 crustal layers
    num = 0
    do ispec = 1, nspec
      if (mesh_data%idoubling(ispec) == IFLAG_CRUST) then
        id = num - (num/3)*3
        mesh_data%idoubling(ispec) = 10 * IFLAG_CRUST + id
        num = num + 1
      endif
    enddo

    if (allocated(xyz)) then
      deallocate(xyz)
    endif
    allocate(xyz(3,num/3))

    ! use only one crustal layer to fine the center
    num = 0
    do ispec = 1, nspec
      if (mesh_data%idoubling(ispec) == 10*IFLAG_CRUST) then
        iglob = mesh_data%ibool(MIDX,MIDY,MIDZ,ispec)
        xyz(:,num+1) = mesh_data%xyz_glob(:,iglob)
        num = num + 1
      endif
    enddo

    ! output lat/lon of the slice center
    xyz_avg = sum(xyz(:,:), dim=2)/num
    call geographic_ecef2ll_zeroalt(xyz_avg(1), xyz_avg(2), xyz_avg(3), lat, lon)

    write(*,"(I06,5X,2F10.3)") iproc, lat*RADIANS_TO_DEGREES, lon*RADIANS_TO_DEGREES

  end do

end program
