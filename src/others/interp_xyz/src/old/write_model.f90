! write out model files for mesh 2

write(*,*) '- write out model files ...'

write(basename,'(a,a,a)') trim(OUTDIR),'/',trim(OUTFNAME)

FMTOUT = "(4E15.5)"

do ipar = 1,NPARS

  fn_model = trim(basename)//'_'//trim(PARNAME(ipar))//'.txt'
  open(unit=IOUT, file=trim(fn_model), status='unknown', iostat=ier, &
       form='formatted')
  if (ier /= 0) then
    print *,'file: ',trim(fn_model)
    stop 'Error opening file'
  end if

  do ipt = 1,npts
    write(IOUT,FMTOUT) xyzi(:,ipt),modeli(ipt,ipar)
  end do

  close(IOUT)

end do