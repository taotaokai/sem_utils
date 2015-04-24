!interpolate mesh data to another mesh
program main

  use constants
  use user_par

  implicit none

  include 'declare_vars.f90'

  ! initialize location data
  uvw = 0.0
  eid = 0
  res = HUGEVAL
  ins = .false.

  ! loop each chunk in mesh 1
  do iproc = 0,nproc-1
    call read_mesh(DATDIR,iproc,NSPEC_MAX,NGLOB_MAX, &
                   nspec,nglob,xyz,xyza,xyzc,ibool,idoubling,model)
    ! locate grid
    include 'locate_gll.f90'
  end do
  
  ! interpolation
  include 'interp_model.f90'
  
  ! write out model file
  include 'write_model.f90'

end program main