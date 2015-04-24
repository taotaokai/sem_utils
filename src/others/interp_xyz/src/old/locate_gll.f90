write(*,*) '- locate grid'
! loop each interpolation points
do ipt = 1,npts

  ! skip point that is already found inside one element of mesh 1
  if (ins(ipt)) cycle 

  ! select elements in source mesh that is close enough to the target point
  idx_sel = (xyzc(1,:)>xyzi(1,ipt)-MARGIN) .and. &
            (xyzc(1,:)<xyzi(1,ipt)+MARGIN) .and. &
            (xyzc(2,:)>xyzi(2,ipt)-MARGIN) .and. &
            (xyzc(2,:)<xyzi(2,ipt)+MARGIN) .and. &
            (xyzc(3,:)>xyzi(3,ipt)-MARGIN) .and. &
            (xyzc(3,:)<xyzi(3,ipt)+MARGIN)

  n_sel = count(idx_sel(1:nspec))
  if (n_sel==0) cycle
  ind_sel(1:n_sel) = pack(IND_MAX,idx_sel)

  ! sort distance
  dist(1:n_sel) = sqrt(sum((xyzc(:,ind_sel(1:n_sel))-xyzi(:,(/(ipt,II=1,n_sel)/)))**2,1))
  call hpsort_index(n_sel,dist,ind_sort)

  ! loop selected elements in mesh 1 to find the one containing the target gll point 
  do i = 1,n_sel
    ispec = ind_sel(ind_sort(i))
    call xyz2cube_bounded(xyza(:,:,ispec),xyzi(:,ipt), uvwi,resi,inside)
    if (inside) then ! record this point and done
      uvw(:,ipt) = uvwi
      eid(ipt) = ispec
      res(ipt) = resi
      ins(ipt) = .true.
      exit
    else ! record the nearest location
      if (resi<res(ipt)) then
        uvw(:,ipt) = uvwi
        eid(ipt) = ispec
        res(ipt) = resi
      ! unset the finish status
        fin(ipt) = .false.
      end if
    end if ! inside
  end do ! i=1,n_sel

end do ! do ipt = 1,npts
