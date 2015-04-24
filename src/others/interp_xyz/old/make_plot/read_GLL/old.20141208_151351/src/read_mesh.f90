! interpolate the GLL model at a given location

subroutine read_mesh &
( iproc, ireg, databasedir, modelname, &
  nspec, nglob, xmesh, ymesh, zmesh, vmesh, ibool)

	implicit none

	include "constants.h"	
	include "values_from_mesher.h"	
	
	!==========================================================
	! specification of subroutine arguments 
	!==========================================================
	
	! inputs
	character(len=*), intent(in) :: databasedir, modelname
	integer, intent(in) :: ireg, iproc
	
	! outputs: mesh data 
	
	! mesh topology 
	integer, intent(out) :: nspec, nglob
	real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE), intent(out) :: xmesh, ymesh, zmesh
	integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), intent(out) :: ibool
	
	! model data
	real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), intent(out) :: vmesh 
	
	!==========================================================
	! define local variables
	!==========================================================
	character(len=256) :: fn_base, fn_topo, fn_model
	
	integer :: ierr ! error value
	
	!==========================================================
	! read topology data
	!==========================================================
	
	! file basename
	write(fn_base,'(a,a,i6.6,a,i1,a)') trim(databasedir),'/proc',iproc,'_reg',ireg,'_'
	
	! read in mesh topology
	fn_topo = trim(fn_base)//'solver_data.bin'
	open(unit=27, file=trim(fn_topo), status='old', action='read', iostat=ierr, form='unformatted')
	if (ierr /= 0) then
	  print*, 'file: ', fn_topo
	  stop 'Error opening file'
	endif
	
	xmesh(:) = 0.0
	ymesh(:) = 0.0
	zmesh(:) = 0.0
	ibool(:,:,:,:) = -1
	read(27) nspec
	read(27) nglob
	read(27) xmesh(1:nglob)
	read(27) ymesh(1:nglob)
	read(27) zmesh(1:nglob)
	read(27) ibool(:,:,:,1:nspec)
	
	close(27)
	
	!==========================================================
	! read model data
	!==========================================================
	fn_model = trim(fn_base)//trim(modelname)//'.bin'
	
	open(unit=27, file=fn_model, status='old', action='read', iostat=ierr, form='unformatted')
	if (ierr /= 0) then
	  print*, 'file: ', fn_model
	  stop 'Error opening file'
	endif
	
	read(27,iostat=ierr) vmesh(:,:,:,1:nspec)
	if (ierr /= 0) then
	  print*, 'file: ', fn_model
	  stop 'Error reading data'
	endif
	
	close(27)
	
end subroutine read_mesh 
