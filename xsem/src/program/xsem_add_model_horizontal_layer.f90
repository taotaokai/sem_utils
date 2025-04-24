subroutine selfdoc()

  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_add_model_horizontal_layer"
  print '(a)', "    - add one horizontal layer into SEM model"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_add_model_LVL1d \"
  print '(a)', "    <mesh_dir> <nproc> <model_dir> <model_tags> "
  print '(a)', "    <layer_r0> <layer_r1> <layer_taperwidth>"
  print '(a)', "    <layer_dlnV> <flag_output_dlnV> <out_dir>"
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', "  The low-velocity layer is a horizontal layer"
  print '(a)', "  the perturbation is tapered near the boundaries"
  print '(a)', "    taper = 0.5 * (1 - cos(PI * d / taper_width)) "
  print '(a)', "  , d is the distance from the layer boundaries"
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (string) mesh_dir:  directory containing proc000***_reg1_solver_data.bin"
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_tag>.bin"
  print '(a)', "  (string) model_tags:  comma delimited string, e.g. vsv,vsh,rho "
  print '(a)', "  (float) layer_r0/r1:  radius range [r0,r1] of plume (km) "
  print '(a)', "  (float) layer_taperwidth:  width of the taper"
  print '(a)', "  (float) layer_dlnV:  velocity perturbation at plume axis"
  print '(a)', "  (int) flag_output_dlnV:  (1/0) 1=output dlnV, 0=don't output dlnV"
  print '(a)', "  (string) out_dir:  output directory of new model"
  print '(a)', ""

end subroutine


!///////////////////////////////////////////////////////////////////////////////
program xsem_add_model_horizontal_layer
  
  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils
  use geographic

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 10
  character(len=MAX_STRING_LEN) :: args(nargs)
  character(len=MAX_STRING_LEN) :: mesh_dir, model_dir, model_tags, out_dir
  integer :: nproc, flag_output_dlnV
  real(dp) :: layer_r0, layer_r1, layer_taperwidth, layer_dlnV

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc

  ! model names
  integer :: nmodel
  character(len=1), parameter :: delimiter = ','
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: iglob, igllx, iglly, igllz, ispec, nspec
  ! model
  real(dp), allocatable :: model_gll(:,:,:,:,:)

  !-- calculate model perturbation
  real(dp) :: xyz(3), r, d
  real(dp) :: taper, layer_r0_add_taper, layer_r1_sub_taper
  real(dp) :: dlnV
  real(dp), allocatable :: dlnV_gll(:,:,:,:)

  !===== read command line arguments

  if (command_argument_count() /= nargs) then
    call selfdoc()
    print *, "[ERROR] xsem_add_model_horizontal_layer: check your input arguments."
    stop
  endif

  do i = 1, nargs
    call get_command_argument(i, args(i))
  enddo
  read(args(1), '(a)') mesh_dir
  read(args(2), *) nproc
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_tags
  read(args(5), *) layer_r0
  read(args(6), *) layer_r1
  read(args(7), *) layer_taperwidth
  read(args(8), *) layer_dlnV
  read(args(9), *) flag_output_dlnV
  read(args(10), '(a)') out_dir

  ! validate input arguments
  if (flag_output_dlnV /= 0 .and. flag_output_dlnV /= 1) then
    print *, "[ERROR] flag_output_dlnV must be 0 or 1"
    stop
  endif

  !===== geometric parameters of plume

  ! change unit: degree -> radians
  ! normalize length by EARTH_R_KM
  layer_r0 = layer_r0 / EARTH_R_KM
  layer_r1 = layer_r1 / EARTH_R_KM
  layer_taperwidth = layer_taperwidth / EARTH_R_KM
  layer_r0_add_taper = layer_r0 + layer_taperwidth
  layer_r1_sub_taper = layer_r1 - layer_taperwidth

  if (layer_r1_sub_taper < layer_r0_add_taper) then
    print *, "[WARN] taper width larger than half layer width"
  endif

  !===== parse model tags

  call sem_utils_delimit_string(model_tags, ',', model_names, nmodel)

  print *, '# nmodel=', nmodel
  print *, '# model_names=', (trim(model_names(i))//"  ", i=1,nmodel) 

  !===== loop each mesh slice
  do iproc = 0, (nproc - 1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)

    nspec = mesh_data%nspec

    ! read model gll
    if (allocated(model_gll)) then
      deallocate(model_gll)
    endif
    allocate(model_gll(nmodel,NGLLX,NGLLY,NGLLZ,nspec))

    if (flag_output_dlnV == 1) then
      if (allocated(dlnV_gll)) then
        deallocate(dlnV_gll)
      endif
      allocate(dlnV_gll(NGLLX,NGLLY,NGLLZ,nspec))
    endif

    call sem_io_read_gll_file_n(model_dir, iproc, iregion, &
                                model_names, nmodel, model_gll)

    ! add slab model on each gll point
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX
    
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)

            ! get the model perturbation
            dlnV = 0.0_dp
            taper = 1.0_dp

            r = sqrt(sum(xyz**2))

            if (r >= layer_r0 .and. r <= layer_r0_add_taper) then
              d = r - layer_r0
              taper = taper * 0.5 * (1 - cos(PI * d / layer_taperwidth))
            endif

            if (r <= layer_r1 .and. r >= layer_r1_sub_taper) then
              d = r - layer_r0
              taper = taper * 0.5 * (1 - cos(PI * d / layer_taperwidth))
            endif

            if (r >= layer_r0 .and. r <= layer_r1) then
              dlnV = taper * layer_dlnV
            endif

            ! get the new model
            model_gll(:,igllx,iglly,igllz,ispec) = (1.0_dp + dlnV) * &
              model_gll(:, igllx,iglly,igllz,ispec)

            if (flag_output_dlnV == 1) then
              dlnV_gll(igllx,iglly,igllz,ispec) = dlnV
            endif

          enddo
        enddo
      enddo
    enddo

    ! write out model
    call sem_io_write_gll_file_n(out_dir, iproc, iregion, &
      model_names, nmodel, model_gll)

    if (flag_output_dlnV == 1) then
      call sem_io_write_gll_file_1(out_dir, iproc, iregion, "dlnV", dlnV_gll)
    endif

  enddo ! iproc

end program
