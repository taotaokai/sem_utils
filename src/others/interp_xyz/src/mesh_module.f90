module mesh

  ! ====================== Specification part ===============

  use constants, only: MAX_STRING_LEN,RK=>CUSTOM_REAL,IIN,IOUT, &
                       HUGEVAL_SNGL, &
                       NDIM,NGLLX,NGLLY,NGLLZ,NGNOD,NGLLCUBE,XGLL, &
                       MESH_SIZE

  implicit none

  ! directories of mesh (_solver_data.bin) and model (_beta_kernel.bin,_vsh.bin, etc.) files
  character(len=MAX_STRING_LEN) :: topo_dir
  character(len=MAX_STRING_LEN) :: model_dir
  integer :: myrank, iregion

  ! mesh geometry
  integer :: nspec, nglob
  real(RK), allocatable :: xyz_glob(:,:)
  real(RK), allocatable :: xyz_anchor(:,:,:)
  real(RK), allocatable :: xyz_center(:,:)
  integer,  allocatable :: ibool(:,:,:,:)

  ! model
  real(RK), allocatable :: model(:,:,:,:,:)
  integer :: num_model
  character(len=MAX_STRING_LEN), allocatable :: model_names(:)

  ! local used variables

  ! locate_xyz
  private :: locate_xyz
  real(RK), allocatable, private :: dist(:)
  logical,  allocatable, private :: idx_sel(:)
  integer,  allocatable, private :: ind_sel(:), ind_sort(:)
  integer,  allocatable, private :: IDX_NSPEC(:)

  ! ====================== Implementation part ===============

contains

! -----------------------
! read mesh geometry
! -----------------------

subroutine read_mesh()

  implicit none

  integer :: ier
  character(len=MAX_STRING_LEN) :: topo_file
  integer, dimension(NGNOD) :: iax, iay, iaz ! index of anchor points
  integer :: ia, ispec, iglob

  write(topo_file,'(a,a,i6.6,a,i1,a)') &
    trim(topo_dir),'/proc',myrank,'_reg',iregion,'_solver_data.bin'

  open(unit=IIN, file=topo_file, status='old', action='read', &
       iostat=ier, form='unformatted')
  if (ier /= 0) then
    print *,'Error open topo file: ', topo_file
    stop
  endif

  ! read mesh dimension
  read(IIN) nspec
  read(IIN) nglob

  ! read in mesh geometry 
  if (.not. allocated(xyz_glob)) then
    allocate(xyz_glob(3,nglob))
  elseif (size(xyz_glob) /= 3*nglob) then
    deallocate(xyz_glob)
    allocate(xyz_glob(3,nglob))
  endif
  if (.not. allocated(ibool)) then
    allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  elseif (size(ibool) /= NGLLCUBE*nspec) then
    deallocate(ibool)
    allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  endif
  read(IIN) xyz_glob(1,1:nglob)
  read(IIN) xyz_glob(2,1:nglob)
  read(IIN) xyz_glob(3,1:nglob)
  read(IIN) ibool(1:NGLLX,1:NGLLY,1:NGLLZ,1:nspec)
  close(IIN)

  ! anchor/center points 
  if (.not. allocated(xyz_anchor)) then
    allocate(xyz_anchor(3,NGNOD,nspec))
  elseif (size(xyz_anchor) /= 3*NGNOD*nspec) then
    deallocate(xyz_anchor)
    allocate(xyz_anchor(3,NGNOD,nspec))
  endif
  if (.not. allocated(xyz_center)) then
    allocate(xyz_center(3,nspec))
  elseif (size(xyz_center) /= 3*nspec) then
    deallocate(xyz_center)
    allocate(xyz_center(3,nspec))
  endif
  call anchor_point_index(iax,iay,iaz)
  ! get anchor/center grids 
  do ispec = 1,nspec
    do ia = 1,NGNOD
      iglob = ibool(iax(ia),iay(ia),iaz(ia),ispec)
      xyz_anchor(:,ia,ispec) = xyz_glob(:,iglob)
    enddo
    ! the last anchor point is the element center
    xyz_center(:,ispec) = xyz_glob(:,iglob)
  enddo

end subroutine


! ----------------------
! locate point (x,y,z)
! ----------------------

subroutine locate_xyz(xyz, uvw, res, eid, iloc)

  implicit none

  real(RK), intent(in) :: xyz(3)
  real(RK), intent(out) :: uvw(3)
  real(RK), intent(out) :: res
  integer, intent(out) :: eid
  integer, intent(out) :: iloc ! 1: inside, 0: close, -1: away of a typical distance

  ! local varaibles
  integer :: i, n_sel, ispec
  real(RK) :: typical_size
  real(RK) :: uvw_i(3), res_i
  logical :: inside

  ! check if mesh data ready
  if ( .not. allocated(xyz_glob) ) then
    stop 'No mesh data read in!'
  endif

  ! select elements in mesh that is close enough to the target point
  !allocate(idx_sel(nspec), ind_sel(nspec), ind_sort(nspec), dist(nspec))
  typical_size = 5*SNGL(MESH_SIZE)
  idx_sel = (xyz_center(1,:)>xyz(1)-typical_size) .and. &
            (xyz_center(1,:)<xyz(1)+typical_size) .and. &
            (xyz_center(2,:)>xyz(2)-typical_size) .and. &
            (xyz_center(2,:)<xyz(2)+typical_size) .and. &
            (xyz_center(3,:)>xyz(3)-typical_size) .and. &
            (xyz_center(3,:)<xyz(3)+typical_size)

  n_sel = count(idx_sel(1:nspec))
  if (n_sel==0) then ! away from the mesh by a distance greater than the element typical_size
    iloc = -1
    return
  endif
  ind_sel(1:n_sel) = pack(IDX_NSPEC,idx_sel)

  ! sort distance
  dist(1:n_sel) = sqrt(sum( (xyz_center(1,ind_sel(1:n_sel)) - xyz(1))**2 + &
                            (xyz_center(2,ind_sel(1:n_sel)) - xyz(2))**2 + &
                            (xyz_center(3,ind_sel(1:n_sel)) - xyz(3))**2 ))
  call hpsort_index(n_sel,dist,ind_sort)

  ! loop selected elements to find the one containing the target xyz 
  res = HUGEVAL_SNGL
  do i = 1, n_sel

    ispec = ind_sel(ind_sort(i))
    call xyz2cube_bounded(xyz_anchor(:,:,ispec),xyz,uvw_i,res_i,inside)
    if (inside) then ! record this point and done
      iloc = 1
      eid = ispec
      uvw = uvw_i
      res = res_i 
      exit
    else ! record the nearest location
      if (res > res_i) then
        iloc = 0
        eid = ispec
        uvw = uvw_i
        res = res_i
      endif
    endif !inside

  enddo

end subroutine


! ----------------------
! read model files
! ----------------------

subroutine read_model()

  implicit none

  integer :: i, ier
  character(len=MAX_STRING_LEN) :: basename, model_file

  if (.not. allocated(model)) then
    allocate(model(NGLLX,NGLLY,NGLLZ,nspec,num_model))
  elseif (size(model) /=  NGLLCUBE*nspec*num_model) then
    deallocate(model)
    allocate(model(NGLLX,NGLLY,NGLLZ,nspec,num_model))
  endif

  write(basename,'(a,a,i6.6,a,i1,a)') &
    trim(model_dir),'/proc',myrank,'_reg',iregion,'_'

  do i = 1, num_model

    model_file = trim(basename)//trim(model_names(i))//'.bin'

    open(unit=IIN, file=trim(model_file), status='old', action='read', &
         iostat=ier, form='unformatted')
    if (ier /= 0) then
      print *, 'Error open model file: ',trim(model_file)
      stop
    end if
    
    read(IIN,iostat=ier) model(:,:,:,1:nspec,i)
    if (ier /= 0) then
      print *, 'Error read model file: ', trim(model_file)
      stop
    endif

    close(IIN)

  enddo

end subroutine


! ----------------------
! interpolate model
! ----------------------

subroutine interp_model(xyz_interp, model_interp, iloc, res_dist)

  implicit none

  real(RK), intent(in) :: xyz_interp(:,:)
  real(RK), intent(inout) :: model_interp(:,:)
  real(RK), intent(inout) :: res_dist(:)
  integer, intent(inout) :: iloc(:)
  
  ! local variables
  integer :: num_interp,i,ix,iy,iz,imodel
  real(RK) :: val, uvw(3), res
  integer :: eid
  real(RK) :: hxi(NGLLX), heta(NGLLY), hgamma(NGLLZ)

  ! prepare store used for locate_xyz
  if (.not. allocated(idx_sel)) then
    allocate(idx_sel(nspec), ind_sel(nspec), ind_sort(nspec), dist(nspec), IDX_NSPEC(nspec))
  elseif (size(idx_sel) /= nspec) then
    deallocate(idx_sel, ind_sel, ind_sort, dist)
    allocate(idx_sel(nspec), ind_sel(nspec), ind_sort(nspec), dist(nspec), IDX_NSPEC(nspec))
  endif
  IDX_NSPEC = (/(i, i=1,nspec)/)

  num_interp = size(xyz_interp,2)
  do i = 1, num_interp

    if (iloc(i) == 1) cycle 

    ! locate point in the current mesh
    call locate_xyz(xyz_interp(:,i), uvw, res, eid, iloc(i))

    if (iloc(i) == 1) then
      res_dist(i) = res
    elseif ( (iloc(i) == 0) .and. (res < res_dist(i))) then
      res_dist(i) = res
    else
      cycle
    endif

    ! interpolate model
    call lagrange_poly(uvw(1),hxi)
    call lagrange_poly(uvw(2),heta)
    call lagrange_poly(uvw(3),hgamma)

    ! loop each parameter
    do imodel = 1, num_model
      val = 0.0
      ! weighted sum of gll value
      do iz = 1, NGLLZ
        do iy = 1, NGLLY
          do ix = 1, NGLLX
            val = val + hxi(ix)*heta(iy)*hgamma(iz)* &
                        model(ix,iy,iz,eid,imodel)
          enddo ! x
        enddo ! y
      enddo ! z
      model_interp(i,imodel) = val
    enddo ! imodel 

  enddo ! i interp

end subroutine


end module
