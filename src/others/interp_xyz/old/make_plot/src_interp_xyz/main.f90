! extract model value at given points xyz
program main

  use constants
  use interp_xyz_par
  
  implicit none

  include 'OUTPUT_FILES/values_from_mesher.f90'

  include 'declare_vars.f90'

  ! read interpolation points from stdin
  include 'read_xyz.f90'

  ! initialize location data
  uvw = 0.0
  eid = 0
  res = HUGEVAL
  ins = .false.
  fin = .false.
	modeli = -1000.0

  ! loop each mesh slice
  do iproc = 0,NPROC-1

    write(*,*) 'iproc=',iproc

    ! read mesh 1
    call read_mesh(DATDIR,iproc,NSPEC_MAX,NGLOB_MAX, &
                   nspec,nglob,xyz,xyza,xyzc,ibool,idoubling,model)

    ! locate grid
    include 'locate_gll.f90'

    ! interpolation
    include 'interp_model.f90'

  enddo ! iproc1
  
  if (any(eid(1:npts)==0)) write(*,*) '- eid=0 detected'
  write(*,*) '- max(res)',maxval(res(1:npts))
 
  ! write out model file
  include 'write_model.f90'

end program main