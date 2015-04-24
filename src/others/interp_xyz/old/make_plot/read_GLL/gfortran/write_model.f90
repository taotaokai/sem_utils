! write out model files for mesh 2

write(*,*) '------ write out model files ...'

write(basename,'(a,a,i6.6,a,i1,a)') trim(OUTDIR2),'/proc',iproc2,'_reg',IREG,'_'

do ipar = 1,NPARS

  fn_model = trim(basename)//trim(PARNAME(ipar))//'.bin'
  open(unit=IOUT, file=trim(fn_model), status='unknown', iostat=ier, form='unformatted')
  if (ier /= 0) then
    print *,'file: ',fn_model
    stop 'Error opening file'
  end if

  write(IOUT) REAL(model2(:,:,:,1:nspec2,ipar))

  close(IOUT)

end do