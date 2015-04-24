! loop each interpolation point
do ipoint = 1,npoint

  ! do nothing if this point is already found inside one mesh element
  if (ins(ipoint)) cycle 

  ! select elements in mesh 1 that is close enough to the target point
  idx_sel = (xyzc(1,:)>xyzi(1,ipoint)-MARGIN) .and. &
            (xyzc(1,:)<xyzi(1,ipoint)+MARGIN) .and. &
            (xyzc(2,:)>xyzi(2,ipoint)-MARGIN) .and. &
            (xyzc(2,:)<xyzi(2,ipoint)+MARGIN) .and. &
            (xyzc(3,:)>xyzi(3,ipoint)-MARGIN) .and. &
            (xyzc(3,:)<xyzi(3,ipoint)+MARGIN)
  n_sel = count(idx_sel(1:nspec))
  if (n_sel==0) cycle
  ind_sel(1:n_sel) = pack(IND_MAX,idx_sel)

  ! sort distance
  dist(1:n_sel) = sqrt(sum((xyzc(:,ind_sel(1:n_sel))-xyzi(:,(/(ipoint,II=1,n_sel)/)))**2,1))
  call hpsort_index(n_sel,dist,ind_sort)

  ! loop selected elements in mesh 1 to find the one containing the target gll point 
  do i = 1,n_sel
    ispec = ind_sel(ind_sort(i))
    call xyz2cube_bounded(xyza(:,:,ispec),xyzi(:,ipoint),uvwi,resi,inside)
    if (inside) then ! record this point and done
      uvw(:,ipoint) = uvwi
      eid(ipoint) = ispec
      res(ipoint) = resi
      ins(ipoint) = .true.
      exit
    else ! record the nearest location
      if (resi<res(ipoint)) then
        uvw(:,ipoint) = uvwi
        eid(ipoint) = ispec
        res(ipoint) = resi
      end if
    end if ! inside
  end do ! i=1,n_sel

end do ! do ipoint = 1,npoint
