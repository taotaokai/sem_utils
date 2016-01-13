subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_combine_solver_data_in_v5.1"
  print '(a)', "    - combine solver_data_1/2 in v5.1 to solver_data in v7.0"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_combine_solver_data_in_v5.1 "
  print '(a)', "    <nproc> <old_mesh_dir> <new_mesh_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) old_mesh_dir:  proc***_reg1_solver_data_1/2.bin v5.1"
  print '(a)', "    and proc***_reg1_array_dims.txt in v5.1"
  print '(a)', "  (string) new_mesh_dir:  output proc***_reg1_solver_data.bin in v7.0"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. can be run in parallel"

end subroutine


program xsem_combine_solver_data

  use sem_constants
  use sem_io

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 3
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: old_mesh_dir
  character(len=MAX_STRING_LEN) :: new_mesh_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc, ier
  ! mesh
  character(len=MAX_STRING_LEN) :: filename
  integer :: nspec, nglob
  real(sp), allocatable :: d1(:), d4(:,:,:,:)
  integer, allocatable :: i1(:), i4(:,:,:,:)
  logical, allocatable :: l1(:)

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    stop "[ERROR] xsem_combine_solver_data: check your inputs."
  endif

  do i = 1, nargs
      call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') old_mesh_dir
  read(args(3), '(a)') new_mesh_dir

  !===== loop each mesh slice
  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! read old mesh array dimensions: reg1_array_dims.txt
    write(filename,"(a,'/proc',i6.6,'_reg',i1,'_array_dims.txt')") &
      trim(old_mesh_dir), iproc, iregion

    open(unit=IIN, file=trim(filename), status='old', iostat=ier)
    if (ier /= 0) then
       write(*,*) '[ERROR] failed to open file: ', trim(filename)
       stop
    endif
    read(IIN,*,iostat=ier) nspec
    if (ier /= 0) stop "[ERROR] failed to read nspec" 
    read(IIN,*,iostat=ier) nglob 
    if (ier /= 0) stop "[ERROR] failed to read nglob" 
    close(IIN)

    print *, 'nspec , nglob = ', nspec, nglob

    ! initialize dummy arrays 
    if (allocated(d1)) then
      deallocate(d1, d4, i1, i4, l1)
    endif
    allocate(d1(nglob), d4(NGLLX,NGLLY,NGLLZ,nspec))
    allocate(i1(nspec), i4(NGLLX,NGLLY,NGLLZ,nspec))
    allocate(l1(nspec))

    ! open new solver_data.bin for write
    call sem_io_open_file_for_write(new_mesh_dir, iproc, iregion, 'solver_data', IOUT)

    ! read old mesh data: solver_data_2.bin
    call sem_io_open_file_for_read(old_mesh_dir, iproc, iregion, 'solver_data_2', IIN)
    ! x coordinate
    read(IIN) d1
    write(IOUT) d1
    ! y coordinate
    read(IIN) d1
    write(IOUT) d1
    ! z coordinate
    read(IIN) d1
    write(IOUT) d1
    ! ibool
    read(IIN) i4
    write(IOUT) i4
    ! idoubling
    read(IIN) i1
    write(IOUT) i1
    ! skip is_on_a_slice_edge (not in v7.0)
    read(IIN) l1
    ! ispec_is_tiso
    read(IIN) l1
    write(IOUT) l1
    ! finish reading solver_data_2.bin
    close(IIN)

    ! read old mesh data: solver_data_1.bin
    call sem_io_open_file_for_read(old_mesh_dir, iproc, iregion, 'solver_data_1', IIN)
    ! xix
    read(IIN) d4
    write(IOUT) d4
    ! xiy
    read(IIN) d4
    write(IOUT) d4
    ! xiz
    read(IIN) d4
    write(IOUT) d4
    ! etax
    read(IIN) d4
    write(IOUT) d4
    ! etay
    read(IIN) d4
    write(IOUT) d4
    ! etaz 
    read(IIN) d4
    write(IOUT) d4
    ! gammax
    read(IIN) d4
    write(IOUT) d4
    ! gammay
    read(IIN) d4
    write(IOUT) d4
    ! gammaz
    read(IIN) d4
    write(IOUT) d4
    ! skip all the rest data because they are not used
    ! finish reading solver_data_1.bin
    close(IIN)

    ! finish writing new mesh solver_data 
    close(IOUT)

  enddo ! iproc

end program
