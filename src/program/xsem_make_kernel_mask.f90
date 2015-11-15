subroutine selfdoc()

    print '(a)', "NAME"
    print '(a)', ""
    print '(a)', "  xsem_make_kernel_mask "
    print '(a)', "    - make kernel mask gll file (mask source, shallow depth)"
    print '(a)', ""
    print '(a)', "SYNOPSIS"
    print '(a)', ""
    print '(a)', "  xsem_make_kernel_mask \"
    print '(a)', "    <nproc> <mesh_dir> <source_vtk_file> <out_dir> \"
    print '(a)', "    <source_mask_radius> <stop_depth> <pass_depth>"
    print '(a)', ""
    print '(a)', "DESCRIPTION"
    print '(a)', ""
    print '(a)', ""
    print '(a)', "PARAMETERS"
    print '(a)', ""
    print '(a)', "  (int) nproc:  number of mesh slices"
    print '(a)', "  (string) mesh_dir:  directory holds proc*_reg1_solver_data.bin"
    print '(a)', "  (string) source_vtk_file:  source.vtk file output from SEM"
    print '(a)', "  (string) out_dir:  output directory for proc*_reg1_source_mask.bin"
    print '(a)', "  (float) source_mask_radius: Gaussian width(one sigma) in km"
    print '(a)', "  (float) stop_depth: stop depth (km) where mask = 0"
    print '(a)', "  (float) pass_depth: pass depth (km) where mask = 1"
    print '(a)', ""

end subroutine


program xsem_make_kernel_mask

    use sem_constants
    use sem_io
    use sem_mesh

    implicit none

    !===== declare variables
    ! command line args
    integer, parameter :: nargs = 7
    character(len=MAX_STRING_LEN) :: args(nargs)
    integer :: nproc
    character(len=MAX_STRING_LEN) :: mesh_dir
    character(len=MAX_STRING_LEN) :: source_vtk_file
    character(len=MAX_STRING_LEN) :: out_dir
    real(dp) :: source_mask_radius
    real(dp) :: stop_depth 
    real(dp) :: pass_depth 

    ! local variables
    integer, parameter :: iregion = IREGION_CRUST_MANTLE ! crust_mantle
    integer :: i, iproc, ier

    ! mesh
    type(sem_mesh_data) :: mesh_data
    integer :: nspec, iglob, igllx, iglly, igllz, ispec

    ! mask gll 
    real(dp), allocatable :: mask(:,:,:,:)
    real(dp) :: xyz(3), weight
    ! source
    character(len=MAX_STRING_LEN) :: dummy 
    real(dp) :: source_xyz(3)
    real(dp) :: dist_sq
    ! depth mask
    real(dp) :: depth

    !===== read command line arguments
    if (command_argument_count() /= nargs) then
        call selfdoc()
        print *, "[ERROR] xsem_make_kernel_mask: check your inputs."
        stop
    endif

    do i = 1, nargs
        call get_command_argument(i, args(i))
    enddo
    read(args(1), *) nproc
    read(args(2), *) mesh_dir
    read(args(3), *) source_vtk_file
    read(args(4), *) out_dir
    read(args(5), *) source_mask_radius
    read(args(6), *) stop_depth
    read(args(7), *) pass_depth
   
    !====== read source.vtk
    open(unit=IIN, file=source_vtk_file)
    ! skip header lines
    read(IIN,*)
    read(IIN,*)
    read(IIN,*)
    read(IIN,*)
    ! nsrc
    read(IIN,*,iostat=ier) dummy, i, dummy
    if (ier /= 0) stop '[ERROR] read source vtk file'
    if (i /= 1) stop '[ERROR] more than one source locations!'
    ! xyz
    read(IIN,*,iostat=ier) source_xyz(1), source_xyz(2), source_xyz(3)
    if (ier /= 0) stop "[ERROR] can't read source location"
    close(IIN)

    print *, "source_xyz: ", source_xyz

    !===== loop each mesh slice
    !non-dimensionalize
    source_mask_radius = source_mask_radius / R_EARTH_KM
    stop_depth = stop_depth / R_EARTH_KM
    pass_depth = pass_depth / R_EARTH_KM

    do iproc = 0, (nproc-1)

        print *, '# iproc=', iproc

        ! read mesh data 
        call sem_mesh_read(mesh_dir, iproc, iregion, mesh_data)
        nspec = mesh_data%nspec

        ! initialize mask gll array 
        if (allocated(mask)) then
            deallocate(mask)
        endif
        allocate(mask(NGLLX,NGLLY,NGLLZ,nspec))

        ! loop each gll point to set mask
        do ispec = 1, nspec
            do igllz = 1, NGLLZ
                do iglly = 1, NGLLY
                    do igllx = 1, NGLLX

                        ! gll point xyz
                        iglob = mesh_data%ibool(igllx,iglly,igllz,ispec)
                        xyz = mesh_data%xyz_glob(:,iglob)
                        depth = 1.0 - sqrt(sum(xyz**2))
                        
                        weight = 1.0_dp
                        ! source mask: mask source region
                        if (source_mask_radius > 0.0) then
                            dist_sq = sum((xyz - source_xyz)**2)
                            weight = weight * (1.0_dp - &
                                exp(-0.5*dist_sq/source_mask_radius**2))
                        endif
                        ! depth mask: mask shallow depth
                        if (pass_depth > stop_depth) then
                            if (stop_depth<depth .and. depth<pass_depth) then
                                weight = weight * (0.5 - 0.5*cos(PI* &
                                    (depth-stop_depth)/(pass_depth-stop_depth)))
                            endif
                        endif

                        mask(igllx,iglly,igllz,ispec) = weight

                    enddo
                enddo
            enddo
        enddo

        print *,'mask value range: ', minval(mask), maxval(mask)

        ! save mask gll file
        call sem_io_write_gll_file_1(out_dir, iproc, iregion, "mask", mask)

    enddo ! iproc

end program
