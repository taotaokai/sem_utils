! read in point coordinates
ipt = 1
do
  read(*,*,iostat=ier) xyzi(:,ipt)
  if (ier/=0) exit
  ipt = ipt+1
  if (ipt>NPTS_MAX) exit
end do
npts = ipt-1

write(*,*) 'npts=',npts
