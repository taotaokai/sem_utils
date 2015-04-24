! check if mesh 1,2 overlap each other

write(*,*) '------ check if two meshes overlap'

idx1_overlap = .true.
idx2_overlap = .true.

n1_overlap = count(idx1_overlap(1:nspec1))
n2_overlap = count(idx2_overlap(1:nspec2))

! loop 
do
  xmin1 = minval(xyzc1(1,:),idx1_overlap); xmax1 = maxval(xyzc1(1,:),idx1_overlap);
  xmin2 = minval(xyzc2(1,:),idx2_overlap); xmax2 = maxval(xyzc2(1,:),idx2_overlap);
  xmin12 = max(xmin1,xmin2); xmax12 = min(xmax1,xmax2);
  
  idx1_overlap = idx1_overlap .and. (xyzc1(1,:)>xmin12-MARGIN) .and. (xyzc1(1,:)<xmax12+MARGIN);
  idx2_overlap = idx2_overlap .and. (xyzc2(1,:)>xmin12-MARGIN) .and. (xyzc2(1,:)<xmax12+MARGIN);

  if ( .not. (any(idx1_overlap).and.any(idx2_overlap))) then
    exit 
  end if
 
  ymin1 = minval(xyzc1(2,:),idx1_overlap); ymax1 = maxval(xyzc1(2,:),idx1_overlap);
  ymin2 = minval(xyzc2(2,:),idx2_overlap); ymax2 = maxval(xyzc2(2,:),idx2_overlap);
  ymin12 = max(ymin1,ymin2); ymax12 = min(ymax1,ymax2);
  
  idx1_overlap = idx1_overlap .and. (xyzc1(2,:)>ymin12-MARGIN) .and. (xyzc1(2,:)<ymax12+MARGIN);
  idx2_overlap = idx2_overlap .and. (xyzc2(2,:)>ymin12-MARGIN) .and. (xyzc2(2,:)<ymax12+MARGIN);

  if ( .not. (any(idx1_overlap).and.any(idx2_overlap))) then
    exit 
  end if
  
  zmin1 = minval(xyzc1(3,:),idx1_overlap); zmax1 = maxval(xyzc1(3,:),idx1_overlap)
  zmin2 = minval(xyzc2(3,:),idx2_overlap); zmax2 = maxval(xyzc2(3,:),idx2_overlap)
  zmin12 = max(zmin1,zmin2); zmax12 = min(zmax1,zmax2)
  
  idx1_overlap = idx1_overlap .and. (xyzc1(3,:)>zmin12-MARGIN) .and. (xyzc1(3,:)<zmax12+MARGIN)
  idx2_overlap = idx2_overlap .and. (xyzc2(3,:)>zmin12-MARGIN) .and. (xyzc2(3,:)<zmax12+MARGIN)

  if ( .not. (any(idx1_overlap).and.any(idx2_overlap))) then
    exit 
  end if

  n1 = count(idx1_overlap(1:nspec1))
  n2 = count(idx2_overlap(1:nspec2))
  
  if ((n1==n1_overlap) .and. (n2==n2_overlap)) then
    exit 
  end if

  n1_overlap = n1
  n2_overlap = n2

end do

if ( .not. (any(idx1_overlap).and.any(idx2_overlap))) then
  write(*,*) '!!!!!! two meshes not overlap!'
  cycle
end if