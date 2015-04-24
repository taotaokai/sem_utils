subroutine get_program_args(databasedir,modelname,nproc,iregion)

	!==========================================================
  ! argument variables
  !==========================================================

	character(len=*) :: databasedir, modelname
	integer :: nproc, iregion

	!==========================================================
  ! local variables
  !==========================================================
	integer :: i
	character(len=256) :: arg(7)

	!==========================================================
  ! read in arguments 
  !==========================================================
  do i = 1, 4
    call get_command_argument(i,arg(i))
    ! usage info
    if (i < 4 .and. len_trim(arg(i)) == 0) then
      print *, ' '
      print *, ' Usage: xinterp_mesh <databasedir> <modelname> <nproc> [<iregion>]'
      stop ' Reenter command line options'
    endif
  enddo

  databasedir = arg(1)
  modelname = arg(2)
	read(arg(3),*) nproc

  ! get region id: 1 for crust_mantle, 2 for outcore, 3 for inner core
  if (len_trim(arg(4)) == 0) then
    iregion  = 1
  else
    read(arg(4),*) iregion
  endif
  if (iregion > 3 .or. iregion <= 0) stop 'Iregion = 1,2,3'

end subroutine get_program_args
