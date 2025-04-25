subroutine selfdoc()
  print '(a)', "NAME"
  print '(a)', ""
  print '(a)', "  xsem_model_gll_reduce_amplitude"
  print '(a)', "    - reduce model amplitude in a given range of depth"
  print '(a)', ""
  print '(a)', "SYNOPSIS"
  print '(a)', ""
  print '(a)', "  xsem_model_gll_reduce_amplitude \"
  print '(a)', "    <nproc> <mesh_dir> <model_dir> <model_name> \"
  print '(a)', "    <depth_stop1> <depth_pass1> <depth_pass2> <depth_stop2> "
  print '(a)', "    <negative_threshold> <positive_threshold> <max_reduce_ratio>"
  print '(a)', "    <out_suffix> <out_dir> "
  print '(a)', ""
  print '(a)', "DESCRIPTION"
  print '(a)', ""
  print '(a)', ""
  print '(a)', "PARAMETERS"
  print '(a)', ""
  print '(a)', "  (int) nproc:  number of mesh slices"
  print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
  print '(a)', "  (string) model_dir:  directory holds proc*_reg1_<model_name>.bin"
  print '(a)', "  (string) model_name:  model name"
  print '(a)', "  (float) depth_stop1/2: depth where the weight decreases to zero"
  print '(a)', "  (float) depth_pass1/2: depth range where the weight is one"
  print '(a)', "  (float) negative_threshold: for model values smaller than this, limit the reduce ratio to <max_reduce_ratio>"
  print '(a)', "  (float) positive_threshold: for model values larger than this, limit the reduce ratio to <max_reduce_ratio>"
  print '(a)', "  (float) max_reduce_ratio: maximum reduction ratio"
  print '(a)', "  (string) out_suffix: proc*_reg1_<model_name><out_suffix>.bin"
  print '(a)', "  (string) out_dir: out_dir/proc*_reg1_<out_tag>.bin"
  print '(a)', ""
  print '(a)', "NOTE"
  print '(a)', ""
  print '(a)', "  1. mostly used for velocity perturbation"
  print '(a)', "  2. weighting_factor = 1.0 - max_reduce_ratio * dlnvs/dlnvs_threshold"

end subroutine


program xsem_model_gll_reduce_amplitude

  use sem_constants
  use sem_io
  use sem_mesh
  use sem_utils

  implicit none

  !===== declare variables

  ! command line args
  integer, parameter :: nargs = 13
  character(len=MAX_STRING_LEN) :: args(nargs)
  integer :: nproc
  character(len=MAX_STRING_LEN) :: mesh_dir
  character(len=MAX_STRING_LEN) :: model_dir
  character(len=MAX_STRING_LEN) :: model_name
  real(dp) :: depth_stop1, depth_pass1
  real(dp) :: depth_stop2, depth_pass2
  real(dp) :: negative_threshold, positive_threshold
  real(dp) :: max_reduce_ratio
  character(len=MAX_STRING_LEN) :: out_suffix
  character(len=MAX_STRING_LEN) :: out_dir

  ! local variables
  integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
  integer :: i, iproc
  ! mesh
  type(sem_mesh_data) :: mesh_data
  integer :: nspec, iglob, igllx, iglly, igllz, ispec
  ! model gll
  real(dp), dimension(:,:,:,:), allocatable :: model_gll
  ! depth taper 
  real(dp) :: xyz(3), taper 
  real(dp) :: depth, depth_width1, depth_width2
  ! reduce factor 
  real(dp) :: model_value, reduce_factor 

  !===== read command line arguments
  if (command_argument_count() /= nargs) then
    call selfdoc()
    stop "[ERROR] check your inputs." 
  endif

  do i = 1, nargs
      call get_command_argument(i, args(i))
  enddo
  read(args(1), *) nproc
  read(args(2), '(a)') mesh_dir
  read(args(3), '(a)') model_dir
  read(args(4), '(a)') model_name 
  read(args(5), *) depth_stop1
  read(args(6), *) depth_pass1
  read(args(7), *) depth_pass2
  read(args(8), *) depth_stop2
  read(args(9), *) negative_threshold 
  read(args(10), *) positive_threshold
  read(args(11), *) max_reduce_ratio 
  read(args(12), '(a)') out_suffix
  read(args(13), '(a)') out_dir

  ! validate inputs
  if ( .not.(depth_stop1<depth_pass1 .and. depth_pass1<depth_pass2 .and. depth_pass2<depth_stop2)) then
    stop "[ERROR] violate condition depth_stop1<depth_pass1<depth_pass2<depth_stop2 "
  endif

  !===== loop each mesh slice

  ! non-dimensionalization
  depth_stop1 = depth_stop1 / EARTH_R_KM
  depth_pass1 = depth_pass1 / EARTH_R_KM
  depth_width1 = depth_pass1 - depth_stop1

  depth_stop2 = depth_stop2 / EARTH_R_KM
  depth_pass2 = depth_pass2 / EARTH_R_KM
  depth_width2 = depth_stop2 - depth_pass2

  do iproc = 0, (nproc-1)

    print *, '# iproc=', iproc

    ! read mesh data 
    call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
    nspec = mesh_data%nspec

    ! intialize arrays 
    if (allocated(model_gll)) then
      deallocate(model_gll)
    end if
    allocate(model_gll(NGLLX,NGLLY,NGLLZ,nspec))

    ! read model 
    call sem_io_read_gll_file_1(model_dir, iproc, iregion, model_name, model_gll)

    ! loop each gll point to set mask
    do ispec = 1, nspec
      do igllz = 1, NGLLZ
        do iglly = 1, NGLLY
          do igllx = 1, NGLLX

            ! gll point xyz
            iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
            xyz = mesh_data%xyz_glob(:,iglob)
            depth = 1.0 - sqrt(sum(xyz**2))

            ! depth taper 
            if (depth <= depth_stop1) then
              taper = 0.0_dp
            elseif (depth >= depth_stop2) then
              taper = 0.0_dp
            elseif (depth_pass1 <= depth .and. depth <= depth_pass2) then
              taper = 1.0_dp
            elseif (depth_stop1 < depth .and. depth < depth_pass1) then
              taper = 0.5 * (1.0 - cos(PI*(depth - depth_stop1)/depth_width1))
            elseif (depth_pass2 < depth .and. depth < depth_stop2) then
              taper = 0.5 * (1.0 + cos(PI*(depth - depth_pass2)/depth_width2))
            endif

            ! reduce factor
            model_value = model_gll(igllx,iglly,igllz,ispec)
            reduce_factor = 0.0
            if (model_value < 0) then
              if (model_value < negative_threshold) then
                model_value = negative_threshold
              endif
              reduce_factor = abs(max_reduce_ratio*model_value/negative_threshold)
            elseif (model_value > 0) then
              if (model_value > positive_threshold) then
                model_value = positive_threshold
              endif
              reduce_factor = abs(max_reduce_ratio*model_value/positive_threshold)
            endif

            ! reduce model amplitude
            model_gll(igllx,iglly,igllz,ispec) = (1.0 - taper*reduce_factor) * model_gll(igllx,iglly,igllz,ispec)

          enddo
        enddo
      enddo
    enddo

    ! save new model gll file
    call sem_io_write_gll_file_1(out_dir, iproc, iregion, trim(model_name)//trim(out_suffix), model_gll)

  enddo ! iproc

end program
